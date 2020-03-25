### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name K562_RNA-seq_10 --genome hg38 --input /project/shefflab/data/guertin/fastq/K562_10pct_RNArc_r2.fastq.gz --single-or-paired single --protocol PRO --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/dev4/results_pipeline -P 4 -M 16000`
*         Compute host:  udc-ba26-28c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/
*  Pipeline started at:   (02-27 09:23:02) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data/guertin/fastq/K562_10pct_RNArc_r2.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_10`
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

Local input file: /project/shefflab/data/guertin/fastq/K562_10pct_RNArc_r2.fastq.gz

> `File_mb`	806.57	PEPPRO	_RES_

> `Read_type`	single	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-27 09:23:02) elapsed: 0.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/K562_10pct_RNArc_r2.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz` (222021)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 222021;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1.fastq`  

> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1.fastq` (222022)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 0.003GB.  
  PID: 222022;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	10000000	PEPPRO	_RES_

> `Fastq_reads`	10000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz']

### FASTQ processing:  (02-27 09:23:31) elapsed: 29.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_processed.fastq`  

> `(cutadapt -j 4 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt` (222072)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 0.169GB.  
  PID: 222072;	Command: cutadapt;	Return code: 0;	Memory used: 0.169GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_processed.fastq` (222110,222111)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 0.169GB.  
  PID: 222110;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 222111;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	8011909.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Uninformative_adapter_reads`	224478.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	2.2448	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/fastqc/K562_RNA-seq_10_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (222143)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.169GB.  
  PID: 222143;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	9775522	PEPPRO	_RES_

> `Trim_loss_rate`	2.24	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_processed.fastq` (222149)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of K562_RNA-seq_10_R1_processed.fastq
Approx 5% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 10% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 15% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 20% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 25% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 30% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 35% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 40% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 45% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 50% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 55% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 60% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 65% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 70% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 75% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 80% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 85% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 90% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 95% complete for K562_RNA-seq_10_R1_processed.fastq
Analysis complete for K562_RNA-seq_10_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 0.187GB.  
  PID: 222149;	Command: fastqc;	Return code: 0;	Memory used: 0.187GB

> `FastQC report r1`	fastqc/K562_RNA-seq_10_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/fastq/processed_R1.flag` (222206)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.187GB.  
  PID: 222206;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Plot adapter insertion distribution (02-27 09:24:57) elapsed: 86.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/cutadapt` (222207)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.221GB.  
  PID: 222207;	Command: Rscript;	Return code: 0;	Memory used: 0.221GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-27 09:25:02) elapsed: 5.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2306	PEPPRO	_RES_

### Prealignments (02-27 09:25:02) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-27 09:25:02) elapsed: 0.0 _TIME_


> `(bowtie2 -p 4 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_10 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
9775522 reads; of these:
  9775522 (100.00%) were unpaired; of these:
    8941856 (91.47%) aligned 0 times
    833666 (8.53%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
8.53% overall alignment rate

> `Aligned_reads_human_rDNA`	833666.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	8.53	PEPPRO	_RES_

### Map to genome (02-27 09:26:15) elapsed: 73.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam`  

> `bowtie2 -p 4 --very-sensitive -X 2000 --rg-id K562_RNA-seq_10 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/tmpmevhvihl -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam` (222518,222519,222520)
<pre>
8941856 reads; of these:
  8941856 (100.00%) were unpaired; of these:
    241092 (2.70%) aligned 0 times
    6440504 (72.03%) aligned exactly 1 time
    2260260 (25.28%) aligned >1 times
97.30% overall alignment rate
[bam_sort_core] merging from 2 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:12:15. Running peak memory: 3.535GB.  
  PID: 222518;	Command: bowtie2;	Return code: 0;	Memory used: 3.535GB  
  PID: 222519;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 222520;	Command: samtools;	Return code: 0;	Memory used: 0.865GB


> `samtools view -q 10 -b -@ 4 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_fail_qc.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam` (223761)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 3.535GB.  
  PID: 223761;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	8700764	PEPPRO	_RES_

> `QC_filtered_reads`	1036653	PEPPRO	_RES_

> `Aligned_reads`	7664111	PEPPRO	_RES_

> `Alignment_rate`	78.4	PEPPRO	_RES_

> `Total_efficiency`	76.64	PEPPRO	_RES_

> `Read_depth`	1.98	PEPPRO	_RES_

### Compress all unmapped read files (02-27 09:46:40) elapsed: 1225.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq` (224787)
<pre>
</pre>
Command completed. Elapsed time: 0:00:50. Running peak memory: 3.535GB.  
  PID: 224787;	Command: pigz;	Return code: 0;	Memory used: 0.004GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam` (224845)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.535GB.  
  PID: 224845;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	227220	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam` (224858)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.535GB.  
  PID: 224858;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_noMT.bam` (224866,224867,224868,224869)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.535GB.  
  PID: 224867;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 224866;	Command: samtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 224868;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 224869;	Command: xargs;	Return code: 0;	Memory used: 0.028GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_noMT.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam` (224897)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 224897;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam` (224898)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.535GB.  
  PID: 224898;	Command: samtools;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-27 09:48:41) elapsed: 121.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -c 4 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv` (224933)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/tmp_K562_RNA-seq_10_sort_a6togo97'
Processing with 4 cores...
Discarding 110 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270731v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 85 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 3.535GB.  
  PID: 224933;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.416GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv`

> `NRF`	0.95	PEPPRO	_RES_

> `PBC1`	0.98	PEPPRO	_RES_

> `PBC2`	50.04	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_unmap.bam`  

> `samtools view -b -@ 4 -f 4  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_unmap.bam` (224965)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.535GB.  
  PID: 224965;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -c -f 4 -@ 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam`

> `Unmapped_reads`	241092	PEPPRO	_RES_

### Split BAM by strand (02-27 09:49:06) elapsed: 24.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam` (224984)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 3.535GB.  
  PID: 224984;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam` (225018)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 3.535GB.  
  PID: 225018;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-27 09:50:24) elapsed: 78.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (225281)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 225281;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/plus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_plus_TssEnrichment.txt` (225283)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.535GB.  
  PID: 225283;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.197GB


> `TSS_coding_score`	15.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/minus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_minus_TssEnrichment.txt` (225300)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.535GB.  
  PID: 225300;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.201GB


> `TSS_non-coding_score`	5.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_minus_TssEnrichment.txt` (225316)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.535GB.  
  PID: 225316;	Command: Rscript;	Return code: 0;	Memory used: 0.292GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_10_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_10_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt` (225333,225334,225335,225336)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 225333;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 225335;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 225334;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 225336;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt` (225338)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 225338;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-27 09:50:43) elapsed: 19.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_tss.bed` (225340,225341)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.535GB.  
  PID: 225340;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 225341;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_gene_body.bed` (225345,225346)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 225346;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 225345;	Command: grep;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSS_density.bed` (225349,225350,225351,225352)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.535GB.  
  PID: 225352;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 225349;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 225351;	Command: sort;	Return code: 0;	Memory used: 0.007GB  
  PID: 225350;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_gene_body_density.bed` (225368,225369,225370)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.535GB.  
  PID: 225368;	Command: bedtools;	Return code: 0;	Memory used: 0.042GB  
  PID: 225370;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 225369;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSS_density.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed` (225384,225385,225386)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 225384;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 225386;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 225385;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	13.88	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed` (225392)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.535GB.  
  PID: 225392;	Command: Rscript;	Return code: 0;	Memory used: 0.33GB

> `Pause index`	QC_hg38/K562_RNA-seq_10_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_10_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed` (225409)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 225409;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (02-27 09:51:20) elapsed: 37.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam`
7664111 2763172

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam`
7664111 2642598

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_gene_sort.bed` (225437,225438)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 225437;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 225438;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_gene_coverage.bed` (225441)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.535GB.  
  PID: 225441;	Command: bedtools;	Return code: 0;	Memory used: 0.041GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed.gz` (225454)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 225454;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed` (225455)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 225455;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (02-27 09:51:55) elapsed: 34.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed` (225473)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.535GB.  
  PID: 225473;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR"` (225476)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 225476;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR_sort.bed` (225477,225478,225479,225480)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 225477;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 225478;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 225480;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 225479;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_plus_coverage.bed` (225482)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 225482;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_minus_coverage.bed` (225489)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 225489;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR"` (225496)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 225496;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR_sort.bed` (225497,225498,225499,225500)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 225497;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 225498;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 225500;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 225499;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_plus_coverage.bed` (225502)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 225502;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_minus_coverage.bed` (225510)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 225510;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer_sort.bed` (225517,225518,225519,225520)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 225517;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 225518;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 225520;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 225519;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_plus_coverage.bed` (225523)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 225523;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_minus_coverage.bed` (225530)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 225530;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon_sort.bed` (225537,225538,225539,225540)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.535GB.  
  PID: 225537;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 225538;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 225540;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 225539;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_plus_coverage.bed` (225544)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.535GB.  
  PID: 225544;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_minus_coverage.bed` (225552)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.535GB.  
  PID: 225552;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron_sort.bed` (225562,225563,225564,225565)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 225562;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 225564;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 225563;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 225565;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_plus_coverage.bed` (225568)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.535GB.  
  PID: 225568;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_minus_coverage.bed` (225577)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.535GB.  
  PID: 225577;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_sort.bed` (225585,225586,225587,225588)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 225585;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 225586;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 225588;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 225587;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_plus_coverage.bed` (225590)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 225590;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_minus_coverage.bed` (225597)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 225597;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region"` (225604)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 225604;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed` (225605,225606,225607,225608)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 225605;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 225607;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 225606;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 225608;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_plus_coverage.bed` (225611)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 225611;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_minus_coverage.bed` (225619)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 225619;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB


### Plot cFRiF/FRiF (02-27 09:53:38) elapsed: 103.0 _TIME_


> `samtools view -@ 4 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_10 -z 3099922541 -n 3790270 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_plus_coverage.bed` (225633)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 3.535GB.  
  PID: 225633;	Command: Rscript;	Return code: 0;	Memory used: 0.511GB

> `cFRiF`	QC_hg38/K562_RNA-seq_10_cFRiF.pdf	cFRiF	QC_hg38/K562_RNA-seq_10_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_10 -z 3099922541 -n 3790270 -y frif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_plus_coverage.bed` (225672)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 3.535GB.  
  PID: 225672;	Command: Rscript;	Return code: 0;	Memory used: 0.511GB

> `FRiF`	QC_hg38/K562_RNA-seq_10_FRiF.pdf	FRiF	QC_hg38/K562_RNA-seq_10_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-27 09:54:33) elapsed: 55.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_exons_sort.bed` (225703,225704)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 225704;	Command: bedtools;	Return code: 0;	Memory used: 0.08GB  
  PID: 225703;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_introns_sort.bed` (225711,225712,225713)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 225712;	Command: bedtools;	Return code: 0;	Memory used: 0.084GB  
  PID: 225711;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 225713;	Command: bedtools;	Return code: 0;	Memory used: 0.092GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_coverage.bed` (225720)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.535GB.  
  PID: 225720;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_coverage.bed` (225736)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.535GB.  
  PID: 225736;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.664111)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_rpkm.bed` (225948,225949,225950)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 225948;	Command: awk;	Return code: 0;	Memory used: 0.006GB  
  PID: 225950;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 225949;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.664111)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_rpkm.bed` (225955,225956,225957)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 225955;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 225957;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 225956;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_rpkm.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed` (225960,225961,225962)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 225960;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 225962;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 225961;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	2.06	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed --annotate` (225968)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.535GB.  
  PID: 225968;	Command: Rscript;	Return code: 0;	Memory used: 0.301GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_10_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_10_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed` (225985)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 225985;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Produce bigWig files (02-27 09:55:17) elapsed: 43.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam` (225992)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.535GB.  
  PID: 225992;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge` (225996)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_10_plus_cuttrace_8115ggs0'
Processing with 1 cores...
stdin is empty of data
Discarding 121 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270719v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_GL000218v1']
Keeping 74 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270718v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrEBV']
Reduce step (merge files)...
Merging 74 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_exact_body_0-mer.bw'
Merging 74 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:14:59. Running peak memory: 3.535GB.  
  PID: 225996;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.5GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam` (228174)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.535GB.  
  PID: 228174;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge` (228178)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_10_minus_cuttrace__tl7etu7'
Processing with 1 cores...
Discarding 124 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270731v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 71 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 71 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_exact_body_0-mer.bw'
Merging 71 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:14:42. Running peak memory: 3.535GB.  
  PID: 228178;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.886GB

Starting cleanup: 52 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  1:02:06
*  Total elapsed time (all runs):  1:16:17
*         Peak memory (this run):  3.5349 GB
*        Pipeline completed time: 2020-02-27 10:25:08
