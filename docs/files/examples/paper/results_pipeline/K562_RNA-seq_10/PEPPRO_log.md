### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name K562_RNA-seq_10 --genome hg38 --input /project/shefflab/data/guertin/fastq/K562_10pct_RNArc_r2.fastq.gz --single-or-paired single --protocol PRO --max-len -1 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/results_pipeline -P 4 -M 16000`
*         Compute host:  udc-ba26-19
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/
*  Pipeline started at:   (02-18 11:08:54) elapsed: 13.0 _TIME_

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
*              `input`:  `['/project/shefflab/data/guertin/fastq/K562_10pct_RNArc_r2.fastq.gz']`
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

### Merge/link and fastq conversion:  (02-18 11:08:56) elapsed: 2.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/K562_10pct_RNArc_r2.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz` (161864)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.002GB.  
  PID: 161864;	Command: ln;	Return code: 0;	Memory used: 0.002GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1.fastq`  

> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1.fastq` (161865)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 0.003GB.  
  PID: 161865;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	10000000	PEPPRO	_RES_

> `Fastq_reads`	10000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz']

### FASTQ processing:  (02-18 11:09:43) elapsed: 48.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_processed.fastq`  

> `(cutadapt -j 4 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt` (161928)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 0.163GB.  
  PID: 161928;	Command: cutadapt;	Return code: 0;	Memory used: 0.163GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_processed.fastq` (162196,162197)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 0.163GB.  
  PID: 162196;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 162197;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	8011909.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_too_short`	224478.0	PEPPRO	_RES_

> `Pct_reads_too_short`	2.2448	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/fastqc/K562_RNA-seq_10_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (162236)
### Calculate the number of trimmed reads
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.163GB.  
  PID: 162236;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	9775522	PEPPRO	_RES_

> `Trim_loss_rate`	2.24	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_processed.fastq` (162242)
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
Command completed. Elapsed time: 0:00:47. Running peak memory: 0.195GB.  
  PID: 162242;	Command: fastqc;	Return code: 0;	Memory used: 0.195GB

> `FastQC report r1`	fastqc/K562_RNA-seq_10_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_

### Plot adapter insertion distribution (02-18 11:11:33) elapsed: 110.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/cutadapt` (162302)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 0.195GB.  
  PID: 162302;	Command: Rscript;	Return code: 0;	Memory used: 0.095GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-18 11:11:43) elapsed: 10.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2306	PEPPRO	_RES_

### Prealignments (02-18 11:11:44) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-18 11:11:44) elapsed: 0.0 _TIME_


> `(bowtie2 -p 4 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_10 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
9775522 reads; of these:
  9775522 (100.00%) were unpaired; of these:
    8941856 (91.47%) aligned 0 times
    833666 (8.53%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
8.53% overall alignment rate

> `Aligned_reads_human_rDNA`	833666.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	8.53	PEPPRO	_RES_

### Map to genome (02-18 11:13:18) elapsed: 95.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam`  

> `bowtie2 -p 4 --very-sensitive -X 2000 --rg-id K562_RNA-seq_10 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/tmp9xc7grcp -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam` (162440,162441,162442)
<pre>
8941856 reads; of these:
  8941856 (100.00%) were unpaired; of these:
    241092 (2.70%) aligned 0 times
    6440504 (72.03%) aligned exactly 1 time
    2260260 (25.28%) aligned >1 times
97.30% overall alignment rate
[bam_sort_core] merging from 2 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:14:37. Running peak memory: 3.535GB.  
  PID: 162441;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 162440;	Command: bowtie2;	Return code: 0;	Memory used: 3.535GB  
  PID: 162442;	Command: samtools;	Return code: 0;	Memory used: 0.865GB


> `samtools view -q 10 -b -@ 4 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam` (164234)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 3.535GB.  
  PID: 164234;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	8700764	PEPPRO	_RES_

> `QC_filtered_reads`	1036653	PEPPRO	_RES_

> `Aligned_reads`	7664111	PEPPRO	_RES_

> `Alignment_rate`	78.4	PEPPRO	_RES_

> `Total_efficiency`	76.64	PEPPRO	_RES_

> `Read_depth`	1.98	PEPPRO	_RES_

### Compress all unmapped read files (02-18 11:36:54) elapsed: 1416.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq` (165206)
<pre>
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 3.535GB.  
  PID: 165206;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam` (165259)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.535GB.  
  PID: 165259;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	227220	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam` (165272)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.535GB.  
  PID: 165272;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_noMT.bam` (165279,165280,165281,165282)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 3.535GB.  
  PID: 165279;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 165281;	Command: grep;	Return code: 0;	Memory used: 0.001GB  
  PID: 165280;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 165282;	Command: xargs;	Return code: 0;	Memory used: 0.026GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam` (165312)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 165312;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam` (165313)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.535GB.  
  PID: 165313;	Command: samtools;	Return code: 0;	Memory used: 0.008GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-18 11:38:53) elapsed: 119.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -c 4 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv` (165348)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/tmp_K562_RNA-seq_10_sort_3xkmmjgx'
Processing with 4 cores...
Discarding 110 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270731v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 85 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 3.535GB.  
  PID: 165348;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.428GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv`

> `NRF`	0.95	PEPPRO	_RES_

> `PBC1`	0.98	PEPPRO	_RES_

> `PBC2`	50.04	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_unmap.bam`  

> `samtools view -b -@ 4 -f 4  /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_unmap.bam` (165383)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.535GB.  
  PID: 165383;	Command: samtools;	Return code: 0;	Memory used: 0.007GB


> `samtools view -c -f 4 -@ 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam`

> `Unmapped_reads`	241092	PEPPRO	_RES_

### Split BAM by strand (02-18 11:39:21) elapsed: 28.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam` (165402)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 3.535GB.  
  PID: 165402;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam` (165534)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 3.535GB.  
  PID: 165534;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-18 11:40:39) elapsed: 78.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (165696)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 165696;	Command: sed;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/plus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_plus_TssEnrichment.txt` (165697)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.535GB.  
  PID: 165697;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.197GB


> `TSS_Plus_Score`	15.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/minus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_minus_TssEnrichment.txt` (165716)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.535GB.  
  PID: 165716;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.201GB


> `TSS_Minus_Score`	5.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_minus_TssEnrichment.txt` (165732)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 165732;	Command: Rscript;	Return code: 0;	Memory used: 0.205GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_10_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_10_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt` (165750,165751,165752,165753)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 165750;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 165752;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 165751;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 165753;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt` (165755)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 165755;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-18 11:41:00) elapsed: 21.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_tss.bed` (165757,165758)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.535GB.  
  PID: 165757;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 165758;	Command: bedtools;	Return code: 0;	Memory used: 0.098GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_gene_body.bed` (165762,165763)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 165763;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB  
  PID: 165762;	Command: grep;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSS_density.bed` (165769,165770,165771,165772)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.535GB.  
  PID: 165772;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 165769;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 165771;	Command: sort;	Return code: 0;	Memory used: 0.009GB  
  PID: 165770;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_gene_body_density.bed` (165785,165786,165787)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.535GB.  
  PID: 165785;	Command: bedtools;	Return code: 0;	Memory used: 0.04GB  
  PID: 165787;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 165786;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed` (165802,165803,165804)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 165802;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 165804;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 165803;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	13.88	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed` (165810)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.535GB.  
  PID: 165810;	Command: Rscript;	Return code: 0;	Memory used: 0.237GB

> `Pause index`	QC_hg38/K562_RNA-seq_10_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_10_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed` (165828)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 165828;	Command: pigz;	Return code: 0;	Memory used: 0.001GB


### Calculate Fraction of Reads in pre-mature mRNA (02-18 11:41:40) elapsed: 39.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam`
7664111 2763172

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam`
7664111 2642598

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_gene_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_gene_sort.bed` (165858,165859)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 165859;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 165858;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_gene_coverage.bed` (165862)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.535GB.  
  PID: 165862;	Command: bedtools;	Return code: 0;	Memory used: 0.041GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed.gz` (165880)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.535GB.  
  PID: 165880;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed` (165884)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 165884;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate fraction and proportion of reads in features (FRiF/PRiF) (02-18 11:42:15) elapsed: 36.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed` (165893)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 165893;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR"` (165895)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 165895;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR_sort.bed` (165896,165897,165898,165899)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 165896;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 165897;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 165899;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 165898;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_plus_coverage.bed` (165901)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 165901;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_minus_coverage.bed` (165909)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 165909;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR"` (165916)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 165916;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR_sort.bed` (165917,165918,165919,165920)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 165917;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 165918;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 165920;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 165919;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_plus_coverage.bed` (165922)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 165922;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_minus_coverage.bed` (165929)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 165929;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer_sort.bed` (165936,165937,165938,165939)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 165936;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 165937;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 165939;	Command: bedtools;	Return code: 0;	Memory used: 0.047GB  
  PID: 165938;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_plus_coverage.bed` (165942)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 165942;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_minus_coverage.bed` (165949)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 165949;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon_sort.bed` (165956,165957,165958,165959)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.535GB.  
  PID: 165956;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 165957;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 165959;	Command: bedtools;	Return code: 0;	Memory used: 0.174GB  
  PID: 165958;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_plus_coverage.bed` (165963)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.535GB.  
  PID: 165963;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_minus_coverage.bed` (165971)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.535GB.  
  PID: 165971;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron_sort.bed` (165983,165984,165985,165986)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.535GB.  
  PID: 165983;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 165985;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 165984;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 165986;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_plus_coverage.bed` (165989)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.535GB.  
  PID: 165989;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_minus_coverage.bed` (165998)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.535GB.  
  PID: 165998;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_sort.bed` (166006,166007,166008,166009)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 166006;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 166007;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 166009;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 166008;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_plus_coverage.bed` (166012)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.535GB.  
  PID: 166012;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_minus_coverage.bed` (166019)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 166019;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region"` (166026)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 166026;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed` (166028,166029,166030,166031)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.535GB.  
  PID: 166028;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 166030;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 166029;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 166031;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_plus_coverage.bed` (166034)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.535GB.  
  PID: 166034;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_minus_coverage.bed` (166041)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.535GB.  
  PID: 166041;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB


### Plot FRiF/PRiF (02-18 11:44:01) elapsed: 106.0 _TIME_


> `samtools view -@ 4 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_10 -z 3099922541 -n 3790270 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_frif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_plus_coverage.bed` (166055)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 3.535GB.  
  PID: 166055;	Command: Rscript;	Return code: 0;	Memory used: 0.528GB

> `FRiF`	QC_hg38/K562_RNA-seq_10_frif.pdf	FRiF	QC_hg38/K562_RNA-seq_10_frif.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_prif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_10 -z 3099922541 -n 3790270 -y prif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_prif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_plus_coverage.bed` (166098)
<pre>
Cumulative prif plot completed!

</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 3.535GB.  
  PID: 166098;	Command: Rscript;	Return code: 0;	Memory used: 0.514GB

> `PRiF`	QC_hg38/K562_RNA-seq_10_prif.pdf	PRiF	QC_hg38/K562_RNA-seq_10_prif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-18 11:45:01) elapsed: 59.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_exons_sort.bed` (166137,166138)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.535GB.  
  PID: 166138;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB  
  PID: 166137;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_introns_sort.bed` (166325,166326,166327)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.535GB.  
  PID: 166326;	Command: bedtools;	Return code: 0;	Memory used: 0.082GB  
  PID: 166325;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 166327;	Command: bedtools;	Return code: 0;	Memory used: 0.098GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_coverage.bed` (166342)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.535GB.  
  PID: 166342;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_coverage.bed` (166354)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.535GB.  
  PID: 166354;	Command: bedtools;	Return code: 0;	Memory used: 0.027GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.664111)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_rpkm.bed` (166368,166369,166370)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 166368;	Command: awk;	Return code: 0;	Memory used: 0.006GB  
  PID: 166370;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 166369;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.664111)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_rpkm.bed` (166372,166373,166374)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.535GB.  
  PID: 166372;	Command: awk;	Return code: 0;	Memory used: 0.006GB  
  PID: 166374;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 166373;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed` (166377,166378,166379)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 166377;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 166379;	Command: sort;	Return code: 0;	Memory used: 0.007GB  
  PID: 166378;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	2.06	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed --annotate` (166385)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.535GB.  
  PID: 166385;	Command: Rscript;	Return code: 0;	Memory used: 0.232GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_10_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_10_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed` (166402)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.535GB.  
  PID: 166402;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Produce bigWig files (02-18 11:45:48) elapsed: 47.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam` (166409)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.535GB.  
  PID: 166409;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge` (166414)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_10_plus_cuttrace__bvyu59m'
Processing with 1 cores...
stdin is empty of data
Discarding 121 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270719v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_GL000218v1']
Keeping 74 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270718v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrEBV']
Reduce step (merge files)...
Merging 74 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_exact_body_0-mer.bw'
Merging 74 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:16:42. Running peak memory: 3.535GB.  
  PID: 166414;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.432GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam` (169294)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.535GB.  
  PID: 169294;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_exact_body_0-mer.bw -p 2 --variable-step --tail-edge` (169329)
<pre>
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_10_minus_cuttrace_2mtnokst'
Processing with 2 cores...
Discarding 124 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270731v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 71 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 71 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_exact_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.535GB.  
  PID: 169329;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.103GB

Starting cleanup: 52 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:54:15
*  Total elapsed time (all runs):  1:11:28
*         Peak memory (this run):  3.535 GB
*        Pipeline completed time: 2020-02-18 12:02:56
