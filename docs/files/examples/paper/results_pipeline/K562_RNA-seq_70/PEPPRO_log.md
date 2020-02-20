### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name K562_RNA-seq_70 --genome hg38 --input /project/shefflab/data/guertin/fastq/K562_70pct_RNArc_r2.fastq.gz --single-or-paired single --protocol PRO --max-len -1 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/results_pipeline -P 4 -M 16000`
*         Compute host:  udc-aj37-14c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/
*  Pipeline started at:   (02-18 11:08:59) elapsed: 0.0 _TIME_

### Version log:

*       Python version:  3.6.5
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.8.9
*        Pipeline hash:  d6fefacd05f25e266e8a162fa20e48a6f6c830d2
*      Pipeline branch:  * dev
*        Pipeline date:  2020-02-17 11:24:59 -0500

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
*              `input`:  `['/project/shefflab/data/guertin/fastq/K562_70pct_RNArc_r2.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `16000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed/peppro/paper/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `False`
*        `sample_name`:  `K562_RNA-seq_70`
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

Local input file: /project/shefflab/data/guertin/fastq/K562_70pct_RNArc_r2.fastq.gz

> `File_mb`	788.96	PEPPRO	_RES_

> `Read_type`	single	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-18 11:09:00) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/raw/K562_RNA-seq_70.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/K562_70pct_RNArc_r2.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/raw/K562_RNA-seq_70.fastq.gz` (83746)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.002GB.  
  PID: 83746;	Command: ln;	Return code: 0;	Memory used: 0.002GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/raw/K562_RNA-seq_70.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1.fastq`  

> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/raw/K562_RNA-seq_70.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1.fastq` (83748)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 0.003GB.  
  PID: 83748;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	10000000	PEPPRO	_RES_

> `Fastq_reads`	10000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/raw/K562_RNA-seq_70.fastq.gz']

### FASTQ processing:  (02-18 11:09:51) elapsed: 51.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1_processed.fastq`  

> `(cutadapt -j 4 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt` (83826)
<pre>
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 0.166GB.  
  PID: 83826;	Command: cutadapt;	Return code: 0;	Memory used: 0.166GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1_processed.fastq` (84118,84119)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 0.166GB.  
  PID: 84118;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 84119;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	5020699.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_too_short`	74592.0	PEPPRO	_RES_

> `Pct_reads_too_short`	0.7459	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/fastqc/K562_RNA-seq_70_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (84161)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.166GB.  
  PID: 84161;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	9925408	PEPPRO	_RES_

> `Trim_loss_rate`	0.75	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1_processed.fastq` (84181)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of K562_RNA-seq_70_R1_processed.fastq
Approx 5% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 10% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 15% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 20% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 25% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 30% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 35% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 40% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 45% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 50% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 55% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 60% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 65% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 70% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 75% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 80% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 85% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 90% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 95% complete for K562_RNA-seq_70_R1_processed.fastq
Analysis complete for K562_RNA-seq_70_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:47. Running peak memory: 0.189GB.  
  PID: 84181;	Command: fastqc;	Return code: 0;	Memory used: 0.189GB

> `FastQC report r1`	fastqc/K562_RNA-seq_70_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_

### Plot adapter insertion distribution (02-18 11:11:58) elapsed: 127.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/cutadapt` (84272)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.205GB.  
  PID: 84272;	Command: Rscript;	Return code: 0;	Memory used: 0.205GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_70_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_70_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-18 11:12:05) elapsed: 7.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2307	PEPPRO	_RES_

### Prealignments (02-18 11:12:05) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-18 11:12:05) elapsed: 0.0 _TIME_


> `(bowtie2 -p 4 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_70 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/prealignments/K562_RNA-seq_70_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
9925408 reads; of these:
  9925408 (100.00%) were unpaired; of these:
    9481694 (95.53%) aligned 0 times
    443714 (4.47%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
4.47% overall alignment rate

> `Aligned_reads_human_rDNA`	443714.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	4.47	PEPPRO	_RES_

### Map to genome (02-18 11:13:28) elapsed: 83.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam`  

> `bowtie2 -p 4 --very-sensitive -X 2000 --rg-id K562_RNA-seq_70 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/prealignments/K562_RNA-seq_70_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/tmpvpuozhy1 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_temp.bam` (84437,84438,84439)
<pre>
9481694 reads; of these:
  9481694 (100.00%) were unpaired; of these:
    1078390 (11.37%) aligned 0 times
    5394074 (56.89%) aligned exactly 1 time
    3009230 (31.74%) aligned >1 times
88.63% overall alignment rate
[bam_sort_core] merging from 2 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:12:37. Running peak memory: 3.533GB.  
  PID: 84438;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 84437;	Command: bowtie2;	Return code: 0;	Memory used: 3.533GB  
  PID: 84439;	Command: samtools;	Return code: 0;	Memory used: 0.865GB


> `samtools view -q 10 -b -@ 4 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam` (85736)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 3.533GB.  
  PID: 85736;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	8403304	PEPPRO	_RES_

> `QC_filtered_reads`	1577377	PEPPRO	_RES_

> `Aligned_reads`	6825927	PEPPRO	_RES_

> `Alignment_rate`	68.77	PEPPRO	_RES_

> `Total_efficiency`	68.26	PEPPRO	_RES_

> `Read_depth`	2.86	PEPPRO	_RES_

### Compress all unmapped read files (02-18 11:32:21) elapsed: 1134.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/prealignments/K562_RNA-seq_70_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/prealignments/K562_RNA-seq_70_human_rDNA_unmap.fq` (86355)
<pre>
</pre>
Command completed. Elapsed time: 0:00:53. Running peak memory: 3.533GB.  
  PID: 86355;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_temp.bam` (86410)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.533GB.  
  PID: 86410;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	477155	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam` (86423)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.533GB.  
  PID: 86423;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_noMT.bam` (86430,86431,86432,86433)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.533GB.  
  PID: 86432;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 86430;	Command: samtools;	Return code: 0;	Memory used: 0.003GB  
  PID: 86431;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 86433;	Command: xargs;	Return code: 0;	Memory used: 0.03GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam` (86459)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 86459;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam` (86460)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.533GB.  
  PID: 86460;	Command: samtools;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-18 11:34:20) elapsed: 119.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -c 4 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_bamQC.tsv` (86492)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/tmp_K562_RNA-seq_70_sort_c5d1tyzk'
Processing with 4 cores...
Discarding 119 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrEBV']
Keeping 76 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 3.533GB.  
  PID: 86492;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.377GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_bamQC.tsv`

> `NRF`	0.83	PEPPRO	_RES_

> `PBC1`	0.94	PEPPRO	_RES_

> `PBC2`	26.62	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_unmap.bam`  

> `samtools view -b -@ 4 -f 4  /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_unmap.bam` (86535)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 86535;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools view -c -f 4 -@ 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_temp.bam`

> `Unmapped_reads`	1078390	PEPPRO	_RES_

### Split BAM by strand (02-18 11:34:48) elapsed: 28.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam` (86557)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.533GB.  
  PID: 86557;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam` (86784)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 3.533GB.  
  PID: 86784;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-18 11:35:55) elapsed: 67.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (86816)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 86816;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/plus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_plus_TssEnrichment.txt` (86818)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.533GB.  
  PID: 86818;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.196GB


> `TSS_Plus_Score`	17.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/minus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_minus_TssEnrichment.txt` (86835)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 86835;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.2GB


> `TSS_Minus_Score`	3.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_minus_TssEnrichment.txt` (86851)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.533GB.  
  PID: 86851;	Command: Rscript;	Return code: 0;	Memory used: 0.201GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_70_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_70_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt` (86872,86873,86874,86875)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 86872;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 86874;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 86873;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 86875;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_keep.txt` (86877)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 86877;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (02-18 11:36:16) elapsed: 21.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_ensembl_tss.bed` (86879,86880)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.533GB.  
  PID: 86879;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 86880;	Command: bedtools;	Return code: 0;	Memory used: 0.094GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_ensembl_gene_body.bed` (86884,86885)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 86884;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 86885;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_TSS_density.bed` (86887,86888,86889,86890)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.533GB.  
  PID: 86890;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 86887;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 86889;	Command: sort;	Return code: 0;	Memory used: 0.008GB  
  PID: 86888;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_gene_body_density.bed` (86910,86911,86912)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.533GB.  
  PID: 86910;	Command: bedtools;	Return code: 0;	Memory used: 0.049GB  
  PID: 86912;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 86911;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.bed` (86928,86929,86930)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 86928;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 86930;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 86929;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	12.41	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.bed` (86935)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 86935;	Command: Rscript;	Return code: 0;	Memory used: 0.215GB

> `Pause index`	QC_hg38/K562_RNA-seq_70_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_70_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.bed` (86954)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 86954;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (02-18 11:36:56) elapsed: 40.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam`
6825927 2637171

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam`
6825927 2621059

> `Minus_FRiP`	0.38	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_gene_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_gene_sort.bed` (86984,86985)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 86984;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 86985;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_gene_coverage.bed` (86988)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.533GB.  
  PID: 86988;	Command: bedtools;	Return code: 0;	Memory used: 0.05GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/raw/hg38_annotations.bed.gz` (87012)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 87012;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/raw/hg38_annotations.bed` (87013)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 87013;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate fraction and proportion of reads in features (FRiF/PRiF) (02-18 11:37:27) elapsed: 31.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/raw/hg38_annotations.bed` (87022)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 87022;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/3_UTR"` (87024)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 87024;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/3_UTR_sort.bed` (87026,87027,87028,87029)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 87026;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 87027;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 87029;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 87028;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_3_UTR_plus_coverage.bed` (87031)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 87031;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_3_UTR_minus_coverage.bed` (87050)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 87050;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/5_UTR"` (87057)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 87057;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/5_UTR_sort.bed` (87058,87059,87060,87061)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 87058;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 87059;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 87061;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 87060;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_5_UTR_plus_coverage.bed` (87063)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 87063;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_5_UTR_minus_coverage.bed` (87070)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 87070;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Enhancer_sort.bed` (87076,87077,87078,87079)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 87076;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 87077;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 87079;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 87078;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Enhancer_plus_coverage.bed` (87082)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 87082;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Enhancer_minus_coverage.bed` (87090)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.533GB.  
  PID: 87090;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Exon_sort.bed` (87097,87099,87100,87101)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.533GB.  
  PID: 87097;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 87099;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 87101;	Command: bedtools;	Return code: 0;	Memory used: 0.162GB  
  PID: 87100;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Exon_plus_coverage.bed` (87111)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.533GB.  
  PID: 87111;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Exon_minus_coverage.bed` (87122)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.533GB.  
  PID: 87122;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Intron_sort.bed` (87130,87131,87132,87133)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 87130;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 87132;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 87131;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 87133;	Command: bedtools;	Return code: 0;	Memory used: 0.08GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Intron_plus_coverage.bed` (87136)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.533GB.  
  PID: 87136;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Intron_minus_coverage.bed` (87156)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.533GB.  
  PID: 87156;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_sort.bed` (87164,87165,87166,87167)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 87164;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 87165;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 87167;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 87166;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_plus_coverage.bed` (87169)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 87169;	Command: bedtools;	Return code: 0;	Memory used: 0.029GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_minus_coverage.bed` (87187)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 87187;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_Flanking_Region"` (87194)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 87194;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed` (87195,87196,87197,87198)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 87195;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 87197;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 87196;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 87198;	Command: bedtools;	Return code: 0;	Memory used: 0.05GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_Flanking_Region_plus_coverage.bed` (87207)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 87207;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_Flanking_Region_minus_coverage.bed` (87213)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 87213;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB


### Plot FRiF/PRiF (02-18 11:39:04) elapsed: 97.0 _TIME_


> `samtools view -@ 4 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_70 -z 3099922541 -n 3233843 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_frif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_Flanking_Region_plus_coverage.bed` (87229)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 3.533GB.  
  PID: 87229;	Command: Rscript;	Return code: 0;	Memory used: 0.486GB

> `FRiF`	QC_hg38/K562_RNA-seq_70_frif.pdf	FRiF	QC_hg38/K562_RNA-seq_70_frif.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_prif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_70 -z 3099922541 -n 3233843 -y prif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_prif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_Flanking_Region_plus_coverage.bed` (87273)
<pre>
Cumulative prif plot completed!

</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.533GB.  
  PID: 87273;	Command: Rscript;	Return code: 0;	Memory used: 0.486GB

> `PRiF`	QC_hg38/K562_RNA-seq_70_prif.pdf	PRiF	QC_hg38/K562_RNA-seq_70_prif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-18 11:40:06) elapsed: 61.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_exons_sort.bed` (87544,87545)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 87545;	Command: bedtools;	Return code: 0;	Memory used: 0.082GB  
  PID: 87544;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_introns_sort.bed` (87557,87558,87559)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 87557;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 87559;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 87558;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exons_coverage.bed` (87578)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.533GB.  
  PID: 87578;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_introns_coverage.bed` (87601)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.533GB.  
  PID: 87601;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/6.825927)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exons_rpkm.bed` (87626,87627,87628)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 87626;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 87628;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 87627;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/6.825927)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_introns_rpkm.bed` (87642,87643,87644)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 87642;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 87644;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 87643;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exon_intron_ratios.bed` (87646,87647,87648)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 87646;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 87648;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 87647;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	6.57	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exon_intron_ratios.bed --annotate` (87655)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.533GB.  
  PID: 87655;	Command: Rscript;	Return code: 0;	Memory used: 0.272GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_70_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_70_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exon_intron_ratios.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exon_intron_ratios.bed` (87679)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 87679;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Produce bigWig files (02-18 11:40:49) elapsed: 43.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam` (87686)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.533GB.  
  PID: 87686;	Command: samtools;	Return code: 0;	Memory used: 0.005GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_plus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge` (87690)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_70_plus_cuttrace_5lflizqy'
Processing with 1 cores...
Discarding 131 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_GL000218v1', 'chrEBV']
Keeping 64 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2']
Reduce step (merge files)...
Merging 64 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_plus_exact_body_0-mer.bw'
Merging 64 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:15:03. Running peak memory: 3.533GB.  
  PID: 87690;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.484GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam` (89902)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.533GB.  
  PID: 89902;	Command: samtools;	Return code: 0;	Memory used: 0.005GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_minus_exact_body_0-mer.bw -p 2 --variable-step --tail-edge` (89906)
<pre>
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_70_minus_cuttrace_c4qoktis'
Processing with 2 cores...
stdin is empty of data
Discarding 129 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr16_KI270728v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 66 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr17_GL000205v2_random', 'chr22_KI270732v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 66 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_minus_exact_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.533GB.  
  PID: 89906;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.125GB

Starting cleanup: 52 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:47:14
*  Total elapsed time (all runs):  1:03:29
*         Peak memory (this run):  3.5329 GB
*        Pipeline completed time: 2020-02-18 11:56:13
