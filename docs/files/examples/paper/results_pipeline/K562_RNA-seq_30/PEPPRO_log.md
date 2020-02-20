### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name K562_RNA-seq_30 --genome hg38 --input /project/shefflab/data/guertin/fastq/K562_30pct_RNArc_r2.fastq.gz --single-or-paired single --protocol PRO --max-len -1 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/results_pipeline -P 4 -M 16000`
*         Compute host:  udc-ba27-8
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/
*  Pipeline started at:   (02-18 11:09:00) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data/guertin/fastq/K562_30pct_RNArc_r2.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_30`
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

Local input file: /project/shefflab/data/guertin/fastq/K562_30pct_RNArc_r2.fastq.gz

> `File_mb`	803.62	PEPPRO	_RES_

> `Read_type`	single	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-18 11:09:00) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/K562_30pct_RNArc_r2.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz` (194119)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.002GB.  
  PID: 194119;	Command: ln;	Return code: 0;	Memory used: 0.002GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1.fastq`  

> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1.fastq` (194120)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 0.003GB.  
  PID: 194120;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	10000000	PEPPRO	_RES_

> `Fastq_reads`	10000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz']

### FASTQ processing:  (02-18 11:09:44) elapsed: 44.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_processed.fastq`  

> `(cutadapt -j 4 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt` (194180)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 0.16GB.  
  PID: 194180;	Command: cutadapt;	Return code: 0;	Memory used: 0.16GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_processed.fastq` (194461,194462)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 0.16GB.  
  PID: 194461;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 194462;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	7013383.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_too_short`	174726.0	PEPPRO	_RES_

> `Pct_reads_too_short`	1.7473	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/fastqc/K562_RNA-seq_30_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (194503)
### Calculate the number of trimmed reads
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.16GB.  
  PID: 194503;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	9825274	PEPPRO	_RES_

> `Trim_loss_rate`	1.75	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_processed.fastq` (194508)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of K562_RNA-seq_30_R1_processed.fastq
Approx 5% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 10% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 15% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 20% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 25% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 30% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 35% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 40% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 45% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 50% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 55% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 60% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 65% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 70% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 75% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 80% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 85% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 90% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 95% complete for K562_RNA-seq_30_R1_processed.fastq
Analysis complete for K562_RNA-seq_30_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:45. Running peak memory: 0.22GB.  
  PID: 194508;	Command: fastqc;	Return code: 0;	Memory used: 0.22GB

> `FastQC report r1`	fastqc/K562_RNA-seq_30_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_

### Plot adapter insertion distribution (02-18 11:11:28) elapsed: 103.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/cutadapt` (194571)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 0.22GB.  
  PID: 194571;	Command: Rscript;	Return code: 0;	Memory used: 0.111GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-18 11:11:37) elapsed: 9.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2303	PEPPRO	_RES_

### Prealignments (02-18 11:11:37) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-18 11:11:37) elapsed: 0.0 _TIME_


> `(bowtie2 -p 4 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_30 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
9825274 reads; of these:
  9825274 (100.00%) were unpaired; of these:
    9122083 (92.84%) aligned 0 times
    703191 (7.16%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
7.16% overall alignment rate

> `Aligned_reads_human_rDNA`	703191.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	7.16	PEPPRO	_RES_

### Map to genome (02-18 11:13:11) elapsed: 94.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam`  

> `bowtie2 -p 4 --very-sensitive -X 2000 --rg-id K562_RNA-seq_30 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/tmpgw9_e5ip -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam` (194881,194882,194883)
<pre>
9122083 reads; of these:
  9122083 (100.00%) were unpaired; of these:
    524877 (5.75%) aligned 0 times
    6088436 (66.74%) aligned exactly 1 time
    2508770 (27.50%) aligned >1 times
94.25% overall alignment rate
[bam_sort_core] merging from 2 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:13:27. Running peak memory: 3.534GB.  
  PID: 194881;	Command: bowtie2;	Return code: 0;	Memory used: 3.534GB  
  PID: 194882;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 194883;	Command: samtools;	Return code: 0;	Memory used: 0.865GB


> `samtools view -q 10 -b -@ 4 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam` (196445)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 3.534GB.  
  PID: 196445;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	8597206	PEPPRO	_RES_

> `QC_filtered_reads`	1214446	PEPPRO	_RES_

> `Aligned_reads`	7382760	PEPPRO	_RES_

> `Alignment_rate`	75.14	PEPPRO	_RES_

> `Total_efficiency`	73.83	PEPPRO	_RES_

> `Read_depth`	2.15	PEPPRO	_RES_

### Compress all unmapped read files (02-18 11:34:56) elapsed: 1305.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq` (856)
<pre>
</pre>
Command completed. Elapsed time: 0:00:48. Running peak memory: 3.534GB.  
  PID: 856;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam` (1119)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.534GB.  
  PID: 1119;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	309258	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam` (1132)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 1132;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_noMT.bam` (1139,1140,1142,1143)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 3.534GB.  
  PID: 1139;	Command: samtools;	Return code: 0;	Memory used: 0.001GB  
  PID: 1142;	Command: grep;	Return code: 0;	Memory used: 0.001GB  
  PID: 1140;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 1143;	Command: xargs;	Return code: 0;	Memory used: 0.027GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam` (1168)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 1168;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam` (1169)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 1169;	Command: samtools;	Return code: 0;	Memory used: 0.008GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-18 11:36:56) elapsed: 120.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -c 4 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv` (1240)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/tmp_K562_RNA-seq_30_sort_ao96ifyw'
Processing with 4 cores...
Discarding 115 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270731v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrEBV']
Keeping 80 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.534GB.  
  PID: 1240;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.379GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv`

> `NRF`	0.93	PEPPRO	_RES_

> `PBC1`	0.97	PEPPRO	_RES_

> `PBC2`	46.16	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_unmap.bam`  

> `samtools view -b -@ 4 -f 4  /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_unmap.bam` (1269)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.534GB.  
  PID: 1269;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools view -c -f 4 -@ 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam`

> `Unmapped_reads`	524877	PEPPRO	_RES_

### Split BAM by strand (02-18 11:37:23) elapsed: 26.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam` (1291)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 3.534GB.  
  PID: 1291;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam` (1334)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 3.534GB.  
  PID: 1334;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-18 11:38:40) elapsed: 77.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (1376)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 1376;	Command: sed;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/plus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_plus_TssEnrichment.txt` (1381)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 1381;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.196GB


> `TSS_Plus_Score`	16.2	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/minus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_minus_TssEnrichment.txt` (1410)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 1410;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.2GB


> `TSS_Minus_Score`	5.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_minus_TssEnrichment.txt` (1426)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 1426;	Command: Rscript;	Return code: 0;	Memory used: 0.223GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_30_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_30_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt` (1444,1445,1446,1447)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 1444;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 1446;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 1445;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 1447;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt` (1449)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 1449;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-18 11:39:00) elapsed: 20.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_tss.bed` (1451,1452)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.534GB.  
  PID: 1451;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 1452;	Command: bedtools;	Return code: 0;	Memory used: 0.098GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_gene_body.bed` (1456,1457)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 1456;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 1457;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSS_density.bed` (1460,1461,1462,1463)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.534GB.  
  PID: 1463;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 1460;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 1462;	Command: sort;	Return code: 0;	Memory used: 0.009GB  
  PID: 1461;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_gene_body_density.bed` (1484,1485,1486)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.534GB.  
  PID: 1484;	Command: bedtools;	Return code: 0;	Memory used: 0.038GB  
  PID: 1486;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 1485;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed` (1504,1505,1506)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 1504;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 1506;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 1505;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	13.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed` (1512)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 1512;	Command: Rscript;	Return code: 0;	Memory used: 0.217GB

> `Pause index`	QC_hg38/K562_RNA-seq_30_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_30_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed` (1542)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 1542;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (02-18 11:39:39) elapsed: 39.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam`
7382760 2719315

> `Plus_FRiP`	0.37	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam`
7382760 2636063

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_gene_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_gene_sort.bed` (1571,1572)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 1572;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 1571;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_gene_coverage.bed` (1575)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.534GB.  
  PID: 1575;	Command: bedtools;	Return code: 0;	Memory used: 0.038GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed.gz` (1842)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 1842;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed` (1843)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 1843;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate fraction and proportion of reads in features (FRiF/PRiF) (02-18 11:40:10) elapsed: 31.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed` (1852)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.534GB.  
  PID: 1852;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR"` (1857)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 1857;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR_sort.bed` (1858,1859,1860,1861)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 1858;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 1859;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 1861;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 1860;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_plus_coverage.bed` (1864)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 1864;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_minus_coverage.bed` (1872)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 1872;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR"` (1879)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 1879;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR_sort.bed` (1880,1881,1882,1883)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 1880;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 1881;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 1883;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 1882;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_plus_coverage.bed` (1886)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 1886;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_minus_coverage.bed` (1893)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 1893;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer_sort.bed` (1903,1904,1905,1906)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 1903;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 1904;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 1906;	Command: bedtools;	Return code: 0;	Memory used: 0.048GB  
  PID: 1905;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_plus_coverage.bed` (1909)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 1909;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_minus_coverage.bed` (1918)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 1918;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon_sort.bed` (1925,1926,1927,1928)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.534GB.  
  PID: 1925;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 1926;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 1928;	Command: bedtools;	Return code: 0;	Memory used: 0.174GB  
  PID: 1927;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_plus_coverage.bed` (1932)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 1932;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_minus_coverage.bed` (1941)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 1941;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron_sort.bed` (1949,1950,1951,1952)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.534GB.  
  PID: 1949;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 1951;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 1950;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 1952;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_plus_coverage.bed` (1955)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.534GB.  
  PID: 1955;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_minus_coverage.bed` (1964)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.534GB.  
  PID: 1964;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_sort.bed` (1985,1986,1987,1988)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 1985;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 1986;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 1988;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 1987;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_plus_coverage.bed` (1994)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 1994;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_minus_coverage.bed` (2012)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2012;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region"` (2043)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2043;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed` (2044,2045,2046,2047)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 2044;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 2046;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 2045;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 2047;	Command: bedtools;	Return code: 0;	Memory used: 0.05GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_plus_coverage.bed` (2053)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2053;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_minus_coverage.bed` (2065)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2065;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB


### Plot FRiF/PRiF (02-18 11:41:59) elapsed: 109.0 _TIME_


> `samtools view -@ 4 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_30 -z 3099922541 -n 3603001 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_frif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_plus_coverage.bed` (2083)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.534GB.  
  PID: 2083;	Command: Rscript;	Return code: 0;	Memory used: 0.53GB

> `FRiF`	QC_hg38/K562_RNA-seq_30_frif.pdf	FRiF	QC_hg38/K562_RNA-seq_30_frif.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_prif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_30 -z 3099922541 -n 3603001 -y prif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_prif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_plus_coverage.bed` (2166)
<pre>
Cumulative prif plot completed!

</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.534GB.  
  PID: 2166;	Command: Rscript;	Return code: 0;	Memory used: 0.529GB

> `PRiF`	QC_hg38/K562_RNA-seq_30_prif.pdf	PRiF	QC_hg38/K562_RNA-seq_30_prif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-18 11:43:02) elapsed: 63.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_exons_sort.bed` (2224,2225)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2225;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB  
  PID: 2224;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_introns_sort.bed` (2234,2235,2236)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2235;	Command: bedtools;	Return code: 0;	Memory used: 0.087GB  
  PID: 2234;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 2236;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_coverage.bed` (2253)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.534GB.  
  PID: 2253;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_coverage.bed` (2275)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.534GB.  
  PID: 2275;	Command: bedtools;	Return code: 0;	Memory used: 0.029GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.38276)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_rpkm.bed` (2300,2301,2302)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 2300;	Command: awk;	Return code: 0;	Memory used: 0.005GB  
  PID: 2302;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 2301;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.38276)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_rpkm.bed` (2305,2306,2307)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 2305;	Command: awk;	Return code: 0;	Memory used: 0.005GB  
  PID: 2307;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 2306;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed` (2309,2310,2311)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2309;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 2311;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 2310;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	3.21	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed --annotate` (2318)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 2318;	Command: Rscript;	Return code: 0;	Memory used: 0.235GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_30_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_30_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed` (2337)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2337;	Command: pigz;	Return code: 0;	Memory used: 0.004GB


### Produce bigWig files (02-18 11:43:51) elapsed: 49.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam` (2344)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.534GB.  
  PID: 2344;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge` (2348)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_30_plus_cuttrace_0dn6lulu'
Processing with 1 cores...
Discarding 126 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_GL000218v1', 'chrEBV']
Keeping 69 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270582v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2']
Reduce step (merge files)...
Merging 69 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_exact_body_0-mer.bw'
Merging 69 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:16:16. Running peak memory: 3.534GB.  
  PID: 2348;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.459GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam` (5096)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.534GB.  
  PID: 5096;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_exact_body_0-mer.bw -p 2 --variable-step --tail-edge` (5100)
<pre>
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_30_minus_cuttrace_0jpicsym'
Processing with 2 cores...
Discarding 125 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr16_KI270728v1_random', 'chr22_KI270731v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 70 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 70 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_exact_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.534GB.  
  PID: 5100;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.109GB

Starting cleanup: 52 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:51:31
*  Total elapsed time (all runs):  1:07:34
*         Peak memory (this run):  3.5344 GB
*        Pipeline completed time: 2020-02-18 12:00:31
