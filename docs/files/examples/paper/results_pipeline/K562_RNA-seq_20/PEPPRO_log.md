### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name K562_RNA-seq_20 --genome hg38 --input /project/shefflab/data/guertin/fastq/K562_20pct_RNArc_r2.fastq.gz --single-or-paired single --protocol PRO --max-len -1 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/results_pipeline -P 4 -M 16000`
*         Compute host:  udc-ba27-8
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/
*  Pipeline started at:   (02-18 11:11:41) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data/guertin/fastq/K562_20pct_RNArc_r2.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_20`
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

Local input file: /project/shefflab/data/guertin/fastq/K562_20pct_RNArc_r2.fastq.gz

> `File_mb`	806.51	PEPPRO	_RES_

> `Read_type`	single	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-18 11:11:41) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/raw/K562_RNA-seq_20.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/K562_20pct_RNArc_r2.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/raw/K562_RNA-seq_20.fastq.gz` (194742)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 194742;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/raw/K562_RNA-seq_20.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/fastq/K562_RNA-seq_20_R1.fastq`  

> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/raw/K562_RNA-seq_20.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/fastq/K562_RNA-seq_20_R1.fastq` (194743)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 0.003GB.  
  PID: 194743;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	10000000	PEPPRO	_RES_

> `Fastq_reads`	10000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/raw/K562_RNA-seq_20.fastq.gz']

### FASTQ processing:  (02-18 11:12:15) elapsed: 34.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/fastq/K562_RNA-seq_20_R1_processed.fastq`  

> `(cutadapt -j 4 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/fastq/K562_RNA-seq_20_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/fastq/K562_RNA-seq_20_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/cutadapt/K562_RNA-seq_20_R1_cutadapt.txt` (194795)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 0.163GB.  
  PID: 194795;	Command: cutadapt;	Return code: 0;	Memory used: 0.163GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/fastq/K562_RNA-seq_20_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/fastq/K562_RNA-seq_20_R1_processed.fastq` (194849,194850)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 0.163GB.  
  PID: 194849;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 194850;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/cutadapt/K562_RNA-seq_20_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	7511497.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/cutadapt/K562_RNA-seq_20_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/cutadapt/K562_RNA-seq_20_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/cutadapt/K562_RNA-seq_20_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_too_short`	199835.0	PEPPRO	_RES_

> `Pct_reads_too_short`	1.9984	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/fastqc/K562_RNA-seq_20_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (194891)
### Calculate the number of trimmed reads
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.163GB.  
  PID: 194891;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	9800165	PEPPRO	_RES_

> `Trim_loss_rate`	2.0	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/fastq/K562_RNA-seq_20_R1_processed.fastq` (194897)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of K562_RNA-seq_20_R1_processed.fastq
Approx 5% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 10% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 15% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 20% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 25% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 30% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 35% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 40% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 45% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 50% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 55% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 60% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 65% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 70% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 75% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 80% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 85% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 90% complete for K562_RNA-seq_20_R1_processed.fastq
Approx 95% complete for K562_RNA-seq_20_R1_processed.fastq
Analysis complete for K562_RNA-seq_20_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 0.193GB.  
  PID: 194897;	Command: fastqc;	Return code: 0;	Memory used: 0.193GB

> `FastQC report r1`	fastqc/K562_RNA-seq_20_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_

### Plot adapter insertion distribution (02-18 11:14:03) elapsed: 108.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/cutadapt/K562_RNA-seq_20_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/cutadapt/K562_RNA-seq_20_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/cutadapt` (194965)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.201GB.  
  PID: 194965;	Command: Rscript;	Return code: 0;	Memory used: 0.201GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_20_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_20_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/cutadapt/K562_RNA-seq_20_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-18 11:14:10) elapsed: 7.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/cutadapt/K562_RNA-seq_20_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/cutadapt/K562_RNA-seq_20_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/cutadapt/K562_RNA-seq_20_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/cutadapt/K562_RNA-seq_20_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/cutadapt/K562_RNA-seq_20_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2305	PEPPRO	_RES_

### Prealignments (02-18 11:14:10) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-18 11:14:10) elapsed: 0.0 _TIME_


> `(bowtie2 -p 4 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_20 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/fastq/K562_RNA-seq_20_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/prealignments/K562_RNA-seq_20_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
9800165 reads; of these:
  9800165 (100.00%) were unpaired; of these:
    9031565 (92.16%) aligned 0 times
    768600 (7.84%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
7.84% overall alignment rate

> `Aligned_reads_human_rDNA`	768600.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	7.84	PEPPRO	_RES_

### Map to genome (02-18 11:15:40) elapsed: 90.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam`  

> `bowtie2 -p 4 --very-sensitive -X 2000 --rg-id K562_RNA-seq_20 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/prealignments/K562_RNA-seq_20_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/tmp9rd16mbs -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_temp.bam` (195309,195310,195311)
<pre>
9031565 reads; of these:
  9031565 (100.00%) were unpaired; of these:
    386898 (4.28%) aligned 0 times
    6261109 (69.32%) aligned exactly 1 time
    2383558 (26.39%) aligned >1 times
95.72% overall alignment rate
[bam_sort_core] merging from 2 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:13:22. Running peak memory: 3.534GB.  
  PID: 195309;	Command: bowtie2;	Return code: 0;	Memory used: 3.534GB  
  PID: 195310;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 195311;	Command: samtools;	Return code: 0;	Memory used: 0.867GB


> `samtools view -q 10 -b -@ 4 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam` (196587)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 3.534GB.  
  PID: 196587;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	8644667	PEPPRO	_RES_

> `QC_filtered_reads`	1123750	PEPPRO	_RES_

> `Aligned_reads`	7520917	PEPPRO	_RES_

> `Alignment_rate`	76.74	PEPPRO	_RES_

> `Total_efficiency`	75.21	PEPPRO	_RES_

> `Read_depth`	2.05	PEPPRO	_RES_

### Compress all unmapped read files (02-18 11:37:39) elapsed: 1319.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/prealignments/K562_RNA-seq_20_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/prealignments/K562_RNA-seq_20_human_rDNA_unmap.fq` (1308)
<pre>
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 3.534GB.  
  PID: 1308;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_temp.bam` (1354)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.534GB.  
  PID: 1354;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	267632	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam` (1370)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 1370;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_noMT.bam` (1394,1395,1396,1397)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.534GB.  
  PID: 1395;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 1394;	Command: samtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 1396;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 1397;	Command: xargs;	Return code: 0;	Memory used: 0.028GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam` (1465)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 1465;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam` (1466)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 1466;	Command: samtools;	Return code: 0;	Memory used: 0.008GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-18 11:39:37) elapsed: 119.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam -c 4 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_bamQC.tsv` (1527)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/tmp_K562_RNA-seq_20_sort_dyto9kp8'
Processing with 4 cores...
Discarding 114 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270731v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 81 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 3.534GB.  
  PID: 1527;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.41GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_bamQC.tsv`

> `NRF`	0.94	PEPPRO	_RES_

> `PBC1`	0.97	PEPPRO	_RES_

> `PBC2`	49.21	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_unmap.bam`  

> `samtools view -b -@ 4 -f 4  /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_unmap.bam` (1583)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.534GB.  
  PID: 1583;	Command: samtools;	Return code: 0;	Memory used: 0.007GB


> `samtools view -c -f 4 -@ 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_temp.bam`

> `Unmapped_reads`	386898	PEPPRO	_RES_

### Split BAM by strand (02-18 11:40:04) elapsed: 27.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_plus.bam` (1834)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 3.534GB.  
  PID: 1834;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_minus.bam` (1914)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 3.534GB.  
  PID: 1914;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-18 11:41:23) elapsed: 79.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (1966)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 1966;	Command: sed;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/plus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_plus_TssEnrichment.txt` (1968)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 1968;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.196GB


> `TSS_Plus_Score`	15.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/minus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_minus_TssEnrichment.txt` (1990)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 1990;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.2GB


> `TSS_Minus_Score`	5.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_minus_TssEnrichment.txt` (2014)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 2014;	Command: Rscript;	Return code: 0;	Memory used: 0.235GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_20_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_20_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt` (2032,2033,2034,2035)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2032;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 2034;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 2033;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 2035;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_keep.txt` (2037)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2037;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-18 11:41:43) elapsed: 20.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_ensembl_tss.bed` (2039,2040)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.534GB.  
  PID: 2039;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 2040;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_ensembl_gene_body.bed` (2049,2050)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2049;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 2050;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_TSS_density.bed` (2055,2056,2057,2058)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.534GB.  
  PID: 2058;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 2055;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 2057;	Command: sort;	Return code: 0;	Memory used: 0.009GB  
  PID: 2056;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_gene_body_density.bed` (2079,2080,2081)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.534GB.  
  PID: 2079;	Command: bedtools;	Return code: 0;	Memory used: 0.04GB  
  PID: 2081;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 2080;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_pause_index.bed` (2108,2109,2110)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2108;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 2110;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 2109;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	13.46	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_pause_index.bed` (2115)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 2115;	Command: Rscript;	Return code: 0;	Memory used: 0.238GB

> `Pause index`	QC_hg38/K562_RNA-seq_20_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_20_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_pause_index.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_pause_index.bed` (2132)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2132;	Command: pigz;	Return code: 0;	Memory used: 0.004GB


### Calculate Fraction of Reads in pre-mature mRNA (02-18 11:42:23) elapsed: 40.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_plus.bam`
7520917 2740027

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_minus.bam`
7520917 2638287

> `Minus_FRiP`	0.35	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_gene_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/signal_hg38/K562_RNA-seq_20_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_gene_sort.bed` (2179,2180)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 2180;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 2179;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/signal_hg38/K562_RNA-seq_20_gene_coverage.bed` (2183)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.534GB.  
  PID: 2183;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/raw/hg38_annotations.bed.gz` (2199)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2199;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/raw/hg38_annotations.bed` (2200)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2200;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate fraction and proportion of reads in features (FRiF/PRiF) (02-18 11:42:55) elapsed: 32.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/raw/hg38_annotations.bed` (2209)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 2209;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/3_UTR"` (2211)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2211;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/3_UTR_sort.bed` (2212,2213,2214,2215)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 2212;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 2213;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 2215;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 2214;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_3_UTR_plus_coverage.bed` (2217)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2217;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_3_UTR_minus_coverage.bed` (2229)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2229;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/5_UTR"` (2241)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2241;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/5_UTR_sort.bed` (2242,2243,2244,2245)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 2242;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 2243;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 2245;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 2244;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_5_UTR_plus_coverage.bed` (2248)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2248;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_5_UTR_minus_coverage.bed` (2257)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2257;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Enhancer_sort.bed` (2264,2265,2266,2267)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 2264;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 2265;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 2267;	Command: bedtools;	Return code: 0;	Memory used: 0.047GB  
  PID: 2266;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Enhancer_plus_coverage.bed` (2270)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2270;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Enhancer_minus_coverage.bed` (2279)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2279;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Exon_sort.bed` (2289,2290,2291,2292)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.534GB.  
  PID: 2289;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 2290;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 2292;	Command: bedtools;	Return code: 0;	Memory used: 0.173GB  
  PID: 2291;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Exon_plus_coverage.bed` (2297)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 2297;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Exon_minus_coverage.bed` (2333)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 2333;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Intron_sort.bed` (2367,2368,2369,2370)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.534GB.  
  PID: 2367;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 2369;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 2368;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 2370;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Intron_plus_coverage.bed` (2373)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.534GB.  
  PID: 2373;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Intron_minus_coverage.bed` (2383)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.534GB.  
  PID: 2383;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter_sort.bed` (2392,2393,2394,2395)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2392;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 2393;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 2395;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 2394;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Promoter_plus_coverage.bed` (2397)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2397;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Promoter_minus_coverage.bed` (2404)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2404;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter_Flanking_Region"` (2412)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2412;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter_Flanking_Region_sort.bed` (2413,2414,2415,2416)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 2413;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 2415;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 2414;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 2416;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Promoter_Flanking_Region_plus_coverage.bed` (2419)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2419;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Promoter_Flanking_Region_minus_coverage.bed` (2429)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2429;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB


### Plot FRiF/PRiF (02-18 11:44:46) elapsed: 110.0 _TIME_


> `samtools view -@ 4 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_20 -z 3099922541 -n 3695852 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_frif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Promoter_Flanking_Region_plus_coverage.bed` (2444)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.534GB.  
  PID: 2444;	Command: Rscript;	Return code: 0;	Memory used: 0.526GB

> `FRiF`	QC_hg38/K562_RNA-seq_20_frif.pdf	FRiF	QC_hg38/K562_RNA-seq_20_frif.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_prif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_20 -z 3099922541 -n 3695852 -y prif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_prif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_Promoter_Flanking_Region_plus_coverage.bed` (2705)
<pre>
Cumulative prif plot completed!

</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.534GB.  
  PID: 2705;	Command: Rscript;	Return code: 0;	Memory used: 0.526GB

> `PRiF`	QC_hg38/K562_RNA-seq_20_prif.pdf	PRiF	QC_hg38/K562_RNA-seq_20_prif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-18 11:45:48) elapsed: 62.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_exons_sort.bed` (2743,2744)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2744;	Command: bedtools;	Return code: 0;	Memory used: 0.084GB  
  PID: 2743;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_introns_sort.bed` (2751,2752,2753)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 2752;	Command: bedtools;	Return code: 0;	Memory used: 0.084GB  
  PID: 2751;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 2753;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_exons_coverage.bed` (2760)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.534GB.  
  PID: 2760;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_introns_coverage.bed` (2773)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.534GB.  
  PID: 2773;	Command: bedtools;	Return code: 0;	Memory used: 0.028GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.520917)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_exons_rpkm.bed` (2796,2797,2798)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 2796;	Command: awk;	Return code: 0;	Memory used: 0.006GB  
  PID: 2798;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 2797;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.520917)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_introns_rpkm.bed` (2800,2801,2802)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 2800;	Command: awk;	Return code: 0;	Memory used: 0.006GB  
  PID: 2802;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 2801;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_exon_intron_ratios.bed` (2805,2806,2807)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2805;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 2807;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 2806;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	2.63	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_exon_intron_ratios.bed --annotate` (2813)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 2813;	Command: Rscript;	Return code: 0;	Memory used: 0.235GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_20_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_20_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_exon_intron_ratios.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/QC_hg38/K562_RNA-seq_20_exon_intron_ratios.bed` (2834)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 2834;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Produce bigWig files (02-18 11:46:36) elapsed: 48.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/signal_hg38/K562_RNA-seq_20_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/signal_hg38/K562_RNA-seq_20_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_plus.bam` (2841)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.534GB.  
  PID: 2841;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/signal_hg38/K562_RNA-seq_20_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/signal_hg38/K562_RNA-seq_20_plus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge` (2846)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_20_plus_cuttrace_8usgb22d'
Processing with 1 cores...
stdin is empty of data
Discarding 124 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_GL000218v1']
Keeping 71 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270582v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrEBV']
Reduce step (merge files)...
Merging 71 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/signal_hg38/K562_RNA-seq_20_plus_exact_body_0-mer.bw'
Merging 71 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/signal_hg38/K562_RNA-seq_20_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:16:08. Running peak memory: 3.534GB.  
  PID: 2846;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.453GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/signal_hg38/K562_RNA-seq_20_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/signal_hg38/K562_RNA-seq_20_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_minus.bam` (6099)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.534GB.  
  PID: 6099;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/signal_hg38/K562_RNA-seq_20_minus_exact_body_0-mer.bw -p 2 --variable-step --tail-edge` (6103)
<pre>
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/aligned_hg38/K562_RNA-seq_20_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_20_minus_cuttrace_6uuqhm40'
Processing with 2 cores...
Discarding 124 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270731v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 71 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 71 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_20/signal_hg38/K562_RNA-seq_20_minus_exact_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.534GB.  
  PID: 6103;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.112GB

Starting cleanup: 52 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:51:31
*  Total elapsed time (all runs):  1:07:17
*         Peak memory (this run):  3.5344 GB
*        Pipeline completed time: 2020-02-18 12:03:11
