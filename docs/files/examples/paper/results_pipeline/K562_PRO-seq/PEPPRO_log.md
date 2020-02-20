### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name K562_PRO-seq --genome hg38 --input /project/shefflab/data/sra_fastq/SRR1554311.fastq.gz /project/shefflab/data/sra_fastq/SRR1554312.fastq.gz --single-or-paired single --protocol PRO --max-len -1 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/results_pipeline -P 16 -M 32000`
*         Compute host:  udc-aj38-13c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/
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
*              `cores`:  `16`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `False`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data/sra_fastq/SRR1554311.fastq.gz', '/project/shefflab/data/sra_fastq/SRR1554312.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `32000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed/peppro/paper/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `False`
*        `sample_name`:  `K562_PRO-seq`
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

Local input file: /project/shefflab/data/sra_fastq/SRR1554311.fastq.gz

> `File_mb`	38517.59	PEPPRO	_RES_

> `Read_type`	single	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-18 11:08:56) elapsed: 2.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz`  

> `cat /project/shefflab/data/sra_fastq/SRR1554311.fastq.gz /project/shefflab/data/sra_fastq/SRR1554312.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz` (116467)
<pre>
</pre>
Command completed. Elapsed time: 0:05:27. Running peak memory: 0.003GB.  
  PID: 116467;	Command: cat;	Return code: 0;	Memory used: 0.003GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1.fastq`  

> `pigz -f -p 16 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1.fastq` (117162)
<pre>
</pre>
Command completed. Elapsed time: 0:19:46. Running peak memory: 0.003GB.  
  PID: 117162;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	496581677	PEPPRO	_RES_

> `Fastq_reads`	496581677	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz']

### FASTQ processing:  (02-18 11:46:53) elapsed: 2278.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_processed.fastq`  

> `(cutadapt -j 16 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt` (121023)
<pre>
</pre>
Command completed. Elapsed time: 0:12:48. Running peak memory: 0.506GB.  
  PID: 121023;	Command: cutadapt;	Return code: 0;	Memory used: 0.506GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_processed.fastq` (122368,122370)
<pre>
</pre>
Command completed. Elapsed time: 0:10:29. Running peak memory: 0.506GB.  
  PID: 122368;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 122370;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	423059183.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_too_short`	12437983.0	PEPPRO	_RES_

> `Pct_reads_too_short`	2.5047	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastqc/K562_PRO-seq_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (124586)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.506GB.  
  PID: 124586;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	484143694	PEPPRO	_RES_

> `Trim_loss_rate`	2.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_processed.fastq` (124804)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of K562_PRO-seq_R1_processed.fastq
Approx 5% complete for K562_PRO-seq_R1_processed.fastq
Approx 10% complete for K562_PRO-seq_R1_processed.fastq
Approx 15% complete for K562_PRO-seq_R1_processed.fastq
Approx 20% complete for K562_PRO-seq_R1_processed.fastq
Approx 25% complete for K562_PRO-seq_R1_processed.fastq
Approx 30% complete for K562_PRO-seq_R1_processed.fastq
Approx 35% complete for K562_PRO-seq_R1_processed.fastq
Approx 40% complete for K562_PRO-seq_R1_processed.fastq
Approx 45% complete for K562_PRO-seq_R1_processed.fastq
Approx 50% complete for K562_PRO-seq_R1_processed.fastq
Approx 55% complete for K562_PRO-seq_R1_processed.fastq
Approx 60% complete for K562_PRO-seq_R1_processed.fastq
Approx 65% complete for K562_PRO-seq_R1_processed.fastq
Approx 70% complete for K562_PRO-seq_R1_processed.fastq
Approx 75% complete for K562_PRO-seq_R1_processed.fastq
Approx 80% complete for K562_PRO-seq_R1_processed.fastq
Approx 85% complete for K562_PRO-seq_R1_processed.fastq
Approx 90% complete for K562_PRO-seq_R1_processed.fastq
Approx 95% complete for K562_PRO-seq_R1_processed.fastq
Analysis complete for K562_PRO-seq_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:29:02. Running peak memory: 0.506GB.  
  PID: 124804;	Command: fastqc;	Return code: 0;	Memory used: 0.305GB

> `FastQC report r1`	fastqc/K562_PRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_

### Plot adapter insertion distribution (02-18 12:47:00) elapsed: 3607.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt` (129032)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.506GB.  
  PID: 129032;	Command: Rscript;	Return code: 0;	Memory used: 0.262GB

> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-18 12:47:04) elapsed: 5.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2312	PEPPRO	_RES_

### Prealignments (02-18 12:47:05) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-18 12:47:05) elapsed: 0.0 _TIME_


> `(bowtie2 -p 16 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_PRO-seq -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/prealignments/K562_PRO-seq_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
484143694 reads; of these:
  484143694 (100.00%) were unpaired; of these:
    439463330 (90.77%) aligned 0 times
    44680364 (9.23%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.23% overall alignment rate

> `Aligned_reads_human_rDNA`	44680364.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.23	PEPPRO	_RES_

### Map to genome (02-18 13:44:21) elapsed: 3437.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam`  

> `bowtie2 -p 16 --very-sensitive -X 2000 --rg-id K562_PRO-seq -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/prealignments/K562_PRO-seq_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/tmpeng2v607 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam` (135376,135382,135383)
<pre>
439463330 reads; of these:
  439463330 (100.00%) were unpaired; of these:
    5393540 (1.23%) aligned 0 times
    328159166 (74.67%) aligned exactly 1 time
    105910624 (24.10%) aligned >1 times
98.77% overall alignment rate
[bam_sort_core] merging from 131 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 3:59:32. Running peak memory: 3.852GB.  
  PID: 135376;	Command: bowtie2;	Return code: 0;	Memory used: 3.852GB  
  PID: 135382;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 135383;	Command: samtools;	Return code: 0;	Memory used: 0.918GB


> `samtools view -q 10 -b -@ 16 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam` (164822)
<pre>
</pre>
Command completed. Elapsed time: 0:11:20. Running peak memory: 3.852GB.  
  PID: 164822;	Command: samtools;	Return code: 0;	Memory used: 0.022GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	434069790	PEPPRO	_RES_

> `QC_filtered_reads`	46741044	PEPPRO	_RES_

> `Aligned_reads`	387328746	PEPPRO	_RES_

> `Alignment_rate`	80.0	PEPPRO	_RES_

> `Total_efficiency`	78.0	PEPPRO	_RES_

> `Read_depth`	29.61	PEPPRO	_RES_

### Compress all unmapped read files (02-18 18:51:39) elapsed: 18438.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/prealignments/K562_PRO-seq_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 16 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/prealignments/K562_PRO-seq_human_rDNA_unmap.fq` (172675)
<pre>
</pre>
Command completed. Elapsed time: 0:11:15. Running peak memory: 3.852GB.  
  PID: 172675;	Command: pigz;	Return code: 0;	Memory used: 0.014GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam` (174287)
<pre>
</pre>
Command completed. Elapsed time: 0:06:31. Running peak memory: 3.852GB.  
  PID: 174287;	Command: samtools;	Return code: 0;	Memory used: 0.027GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	9149071	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam` (177099)
<pre>
</pre>
Command completed. Elapsed time: 0:05:45. Running peak memory: 3.852GB.  
  PID: 177099;	Command: samtools;	Return code: 0;	Memory used: 0.028GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 16 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_noMT.bam` (178029,178030,178031,178032)
<pre>
</pre>
Command completed. Elapsed time: 0:05:44. Running peak memory: 3.852GB.  
  PID: 178031;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 178029;	Command: samtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 178030;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 178032;	Command: xargs;	Return code: 0;	Memory used: 0.161GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam` (178777)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.852GB.  
  PID: 178777;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam` (178785)
<pre>
</pre>
Command completed. Elapsed time: 0:05:49. Running peak memory: 3.852GB.  
  PID: 178785;	Command: samtools;	Return code: 0;	Memory used: 0.027GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-18 19:40:40) elapsed: 2940.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -c 16 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv` (181567)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/tmp_K562_PRO-seq_sort_2jmgy6cx'
Processing with 16 cores...
Discarding 81 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1']
Keeping 114 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270411v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:06:20. Running peak memory: 13.655GB.  
  PID: 181567;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 13.655GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv`

> `NRF`	0.28	PEPPRO	_RES_

> `PBC1`	0.54	PEPPRO	_RES_

> `PBC2`	3.89	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_unmap.bam`  

> `samtools view -b -@ 16 -f 4  /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_unmap.bam` (182230)
<pre>
</pre>
Command completed. Elapsed time: 0:01:16. Running peak memory: 13.655GB.  
  PID: 182230;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools view -c -f 4 -@ 16 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam`

> `Unmapped_reads`	5393540	PEPPRO	_RES_

### Split BAM by strand (02-18 19:49:25) elapsed: 526.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam` (182421)
<pre>
</pre>
Command completed. Elapsed time: 0:29:38. Running peak memory: 13.655GB.  
  PID: 182421;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam` (186204)
<pre>
</pre>
Command completed. Elapsed time: 0:27:03. Running peak memory: 13.655GB.  
  PID: 186204;	Command: samtools;	Return code: 0;	Memory used: 0.007GB


### Calculate TSS enrichment (02-18 20:46:06) elapsed: 3401.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (189787)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 13.655GB.  
  PID: 189787;	Command: sed;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/plus_TSS.tsv -p ends -c 16 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_plus_TssEnrichment.txt` (189793)
<pre>
</pre>
Command completed. Elapsed time: 0:00:47. Running peak memory: 13.655GB.  
  PID: 189793;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.268GB


> `TSS_Plus_Score`	13.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/minus_TSS.tsv -p ends -c 16 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_minus_TssEnrichment.txt` (189886)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 13.655GB.  
  PID: 189886;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.781GB


> `TSS_Minus_Score`	5.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_minus_TssEnrichment.txt` (189980)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 13.655GB.  
  PID: 189980;	Command: Rscript;	Return code: 0;	Memory used: 0.205GB

> `TSS enrichment`	QC_hg38/K562_PRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt` (190007,190008,190009,190010)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 13.655GB.  
  PID: 190007;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 190009;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 190008;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 190010;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt` (190012)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 13.655GB.  
  PID: 190012;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-18 20:47:38) elapsed: 92.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_tss.bed` (190014,190015)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 13.655GB.  
  PID: 190014;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 190015;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed` (190018,190019)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 13.655GB.  
  PID: 190018;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 190019;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_TSS_density.bed` (190022,190023,190024,190025)
<pre>
</pre>
Command completed. Elapsed time: 0:09:55. Running peak memory: 13.655GB.  
  PID: 190025;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 190022;	Command: bedtools;	Return code: 0;	Memory used: 0.103GB  
  PID: 190024;	Command: sort;	Return code: 0;	Memory used: 0.012GB  
  PID: 190023;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_gene_body_density.bed` (191277,191282,191283)
<pre>
</pre>
Command completed. Elapsed time: 0:17:07. Running peak memory: 13.655GB.  
  PID: 191282;	Command: awk;	Return code: 0;	Memory used: 0.002GB  
  PID: 191277;	Command: bedtools;	Return code: 0;	Memory used: 1.019GB  
  PID: 191283;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed` (193149,193155,193156)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 13.655GB.  
  PID: 193149;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 193156;	Command: env;	Return code: 0;	Memory used: 0.007GB  
  PID: 193155;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed` (193162)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 13.655GB.  
  PID: 193162;	Command: Rscript;	Return code: 0;	Memory used: 0.235GB

> `Pause index`	QC_hg38/K562_PRO-seq_pause_index.pdf	Pause index	QC_hg38/K562_PRO-seq_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed.gz`  

> `pigz -f -p 16 -f /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed` (193185)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 13.655GB.  
  PID: 193185;	Command: pigz;	Return code: 0;	Memory used: 0.008GB


### Calculate Fraction of Reads in pre-mature mRNA (02-18 21:14:48) elapsed: 1630.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam`
387328746 138157932

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam`
387328746 131399421

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_gene_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_gene_sort.bed` (194632,194634)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 13.655GB.  
  PID: 194632;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 194634;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_gene_coverage.bed` (194636)
<pre>
</pre>
Command completed. Elapsed time: 0:16:27. Running peak memory: 13.655GB.  
  PID: 194636;	Command: bedtools;	Return code: 0;	Memory used: 1.019GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed.gz` (196864)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 13.655GB.  
  PID: 196864;	Command: ln;	Return code: 0;	Memory used: 0.001GB


> `pigz -f -p 16 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed` (196881)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 13.655GB.  
  PID: 196881;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate fraction and proportion of reads in features (FRiF/PRiF) (02-18 21:42:59) elapsed: 1691.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed` (196890)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 13.655GB.  
  PID: 196890;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR"` (196892)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 13.655GB.  
  PID: 196892;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR_sort.bed` (196893,196894,196895,196896)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 13.655GB.  
  PID: 196893;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 196894;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 196896;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 196895;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_plus_coverage.bed` (196898)
<pre>
</pre>
Command completed. Elapsed time: 0:05:21. Running peak memory: 13.655GB.  
  PID: 196898;	Command: bedtools;	Return code: 0;	Memory used: 0.105GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_minus_coverage.bed` (197463)
<pre>
</pre>
Command completed. Elapsed time: 0:05:07. Running peak memory: 13.655GB.  
  PID: 197463;	Command: bedtools;	Return code: 0;	Memory used: 0.162GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR"` (198159)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 13.655GB.  
  PID: 198159;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR_sort.bed` (198160,198161,198162,198163)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 13.655GB.  
  PID: 198160;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 198161;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 198163;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 198162;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_plus_coverage.bed` (198166)
<pre>
</pre>
Command completed. Elapsed time: 0:05:14. Running peak memory: 13.655GB.  
  PID: 198166;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_minus_coverage.bed` (198710)
<pre>
</pre>
Command completed. Elapsed time: 0:05:00. Running peak memory: 13.655GB.  
  PID: 198710;	Command: bedtools;	Return code: 0;	Memory used: 0.064GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer_sort.bed` (199295,199298,199299,199300)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 13.655GB.  
  PID: 199295;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 199298;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 199300;	Command: bedtools;	Return code: 0;	Memory used: 0.048GB  
  PID: 199299;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_plus_coverage.bed` (199303)
<pre>
</pre>
Command completed. Elapsed time: 0:04:52. Running peak memory: 13.655GB.  
  PID: 199303;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_minus_coverage.bed` (199903)
<pre>
</pre>
Command completed. Elapsed time: 0:04:48. Running peak memory: 13.655GB.  
  PID: 199903;	Command: bedtools;	Return code: 0;	Memory used: 0.069GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Exon_sort.bed` (200417,200418,200419,200420)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 13.655GB.  
  PID: 200417;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 200418;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 200420;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 200419;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_plus_coverage.bed` (200426)
<pre>
</pre>
Command completed. Elapsed time: 0:06:16. Running peak memory: 13.655GB.  
  PID: 200426;	Command: bedtools;	Return code: 0;	Memory used: 0.466GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_minus_coverage.bed` (201067)
<pre>
</pre>
Command completed. Elapsed time: 0:05:44. Running peak memory: 13.655GB.  
  PID: 201067;	Command: bedtools;	Return code: 0;	Memory used: 0.185GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Intron_sort.bed` (201970,201971,201972,201973)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 13.655GB.  
  PID: 201970;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 201972;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 201971;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 201973;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_plus_coverage.bed` (201976)
<pre>
</pre>
Command completed. Elapsed time: 0:07:16. Running peak memory: 13.655GB.  
  PID: 201976;	Command: bedtools;	Return code: 0;	Memory used: 0.45GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_minus_coverage.bed` (202715)
<pre>
</pre>
Command completed. Elapsed time: 0:07:08. Running peak memory: 13.655GB.  
  PID: 202715;	Command: bedtools;	Return code: 0;	Memory used: 0.433GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_sort.bed` (203423,203428,203429,203430)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 13.655GB.  
  PID: 203423;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 203428;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 203430;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 203429;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_plus_coverage.bed` (203433)
<pre>
</pre>
Command completed. Elapsed time: 0:05:33. Running peak memory: 13.655GB.  
  PID: 203433;	Command: bedtools;	Return code: 0;	Memory used: 0.582GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_minus_coverage.bed` (204210)
<pre>
</pre>
Command completed. Elapsed time: 0:05:26. Running peak memory: 13.655GB.  
  PID: 204210;	Command: bedtools;	Return code: 0;	Memory used: 0.474GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region"` (204824)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 13.655GB.  
  PID: 204824;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed` (204830,204831,204832,204833)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 13.655GB.  
  PID: 204830;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 204832;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 204831;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 204833;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (204836)
<pre>
</pre>
Command completed. Elapsed time: 0:05:15. Running peak memory: 13.655GB.  
  PID: 204836;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_minus_coverage.bed` (205386)
<pre>
</pre>
Command completed. Elapsed time: 0:05:08. Running peak memory: 13.655GB.  
  PID: 205386;	Command: bedtools;	Return code: 0;	Memory used: 0.166GB


### Plot FRiF/PRiF (02-18 23:01:20) elapsed: 4702.0 _TIME_


> `samtools view -@ 16 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_PRO-seq -z 3099922541 -n 192750289 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_frif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (206012)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 13.655GB.  
  PID: 206012;	Command: Rscript;	Return code: 0;	Memory used: 0.521GB

> `FRiF`	QC_hg38/K562_PRO-seq_frif.pdf	FRiF	QC_hg38/K562_PRO-seq_frif.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_prif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_PRO-seq -z 3099922541 -n 192750289 -y prif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_prif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (206077)
<pre>
Cumulative prif plot completed!

</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 13.655GB.  
  PID: 206077;	Command: Rscript;	Return code: 0;	Memory used: 0.512GB

> `PRiF`	QC_hg38/K562_PRO-seq_prif.pdf	PRiF	QC_hg38/K562_PRO-seq_prif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-18 23:03:01) elapsed: 101.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_exons_sort.bed` (206122,206123)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 13.655GB.  
  PID: 206123;	Command: bedtools;	Return code: 0;	Memory used: 0.094GB  
  PID: 206122;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_introns_sort.bed` (206134,206135,206136)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 13.655GB.  
  PID: 206134;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 206136;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 206135;	Command: bedtools;	Return code: 0;	Memory used: 0.031GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_coverage.bed` (206143)
<pre>
</pre>
Command completed. Elapsed time: 0:10:54. Running peak memory: 13.655GB.  
  PID: 206143;	Command: bedtools;	Return code: 0;	Memory used: 0.292GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_coverage.bed` (207347)
<pre>
</pre>
Command completed. Elapsed time: 0:15:30. Running peak memory: 13.655GB.  
  PID: 207347;	Command: bedtools;	Return code: 0;	Memory used: 0.486GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/387.328746)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_rpkm.bed` (209079,209084,209085)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 13.655GB.  
  PID: 209079;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 209085;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 209084;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/387.328746)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_rpkm.bed` (209088,209089,209090)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 13.655GB.  
  PID: 209088;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 209090;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 209089;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed` (209093,209094,209095)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 13.655GB.  
  PID: 209093;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 209095;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 209094;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed --annotate` (209101)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 13.655GB.  
  PID: 209101;	Command: Rscript;	Return code: 0;	Memory used: 0.234GB

> `mRNA contamination`	QC_hg38/K562_PRO-seq_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_PRO-seq_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed.gz`  

> `pigz -f -p 16 -f /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed` (209124)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 13.655GB.  
  PID: 209124;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (02-18 23:29:44) elapsed: 1603.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam` (209133)
<pre>
</pre>
Command completed. Elapsed time: 0:02:57. Running peak memory: 13.655GB.  
  PID: 209133;	Command: samtools;	Return code: 0;	Memory used: 0.02GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_smooth_body_0-mer.bw -p 10 --variable-step --tail-edge` (209556)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_plus_cuttrace_ikxaws6d'
Processing with 5 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 90 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1']
Keeping 105 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 105 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_exact_body_0-mer.bw'
Merging 105 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:13:25. Running peak memory: 13.655GB.  
  PID: 209556;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.474GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam` (212099)
<pre>
</pre>
Command completed. Elapsed time: 0:02:49. Running peak memory: 13.655GB.  
  PID: 212099;	Command: samtools;	Return code: 0;	Memory used: 0.021GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_exact_body_0-mer.bw -p 10 --variable-step --tail-edge` (212277)
<pre>
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_minus_cuttrace_259r3xl_'
Processing with 10 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 93 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 102 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270411v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270336v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 102 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_exact_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:03:00. Running peak memory: 13.655GB.  
  PID: 212277;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.656GB

Starting cleanup: 52 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  12:43:25
*  Total elapsed time (all runs):  17:21:05
*         Peak memory (this run):  13.6551 GB
*        Pipeline completed time: 2020-02-18 23:52:06
