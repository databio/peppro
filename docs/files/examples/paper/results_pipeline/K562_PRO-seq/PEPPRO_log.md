### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name K562_PRO-seq --genome hg38 --input /project/shefflab/data/sra_fastq/SRR1554311.fastq.gz /project/shefflab/data/sra_fastq/SRR1554312.fastq.gz --single-or-paired single --protocol PRO --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/dev4/results_pipeline -P 16 -M 32000`
*         Compute host:  udc-aj40-18c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/
*  Pipeline started at:   (02-27 09:23:05) elapsed: 2.0 _TIME_

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
*      `output_parent`:  `/project/shefflab/processed/peppro/paper/dev4/results_pipeline`
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

### Merge/link and fastq conversion:  (02-27 09:23:05) elapsed: 0.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz`  

> `cat /project/shefflab/data/sra_fastq/SRR1554311.fastq.gz /project/shefflab/data/sra_fastq/SRR1554312.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz` (290725)
<pre>
</pre>
Command completed. Elapsed time: 0:03:23. Running peak memory: 0.003GB.  
  PID: 290725;	Command: cat;	Return code: 0;	Memory used: 0.003GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1.fastq`  

> `pigz -f -p 16 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1.fastq` (291139)
<pre>
</pre>
Command completed. Elapsed time: 0:19:04. Running peak memory: 0.003GB.  
  PID: 291139;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	496581677	PEPPRO	_RES_

> `Fastq_reads`	496581677	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz']

### FASTQ processing:  (02-27 09:59:18) elapsed: 2173.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_processed.fastq`  

> `(cutadapt -j 16 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt` (294456)
<pre>
</pre>
Command completed. Elapsed time: 0:14:37. Running peak memory: 0.49GB.  
  PID: 294456;	Command: cutadapt;	Return code: 0;	Memory used: 0.49GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_processed.fastq` (296075,296076)
<pre>
</pre>
Command completed. Elapsed time: 0:09:28. Running peak memory: 0.49GB.  
  PID: 296076;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB  
  PID: 296075;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	423059183.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Uninformative_adapter_reads`	12437983.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	2.5047	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/fastqc/K562_PRO-seq_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (297669)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.49GB.  
  PID: 297669;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	484143694	PEPPRO	_RES_

> `Trim_loss_rate`	2.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_processed.fastq` (298018)
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
Command completed. Elapsed time: 0:27:52. Running peak memory: 0.49GB.  
  PID: 298018;	Command: fastqc;	Return code: 0;	Memory used: 0.308GB

> `FastQC report r1`	fastqc/K562_PRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/fastq/processed_R1.flag` (300791)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.49GB.  
  PID: 300791;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (02-27 10:59:23) elapsed: 3604.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/cutadapt` (300797)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.49GB.  
  PID: 300797;	Command: Rscript;	Return code: 0;	Memory used: 0.202GB

> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-27 10:59:29) elapsed: 7.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2312	PEPPRO	_RES_

### Prealignments (02-27 10:59:29) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-27 10:59:29) elapsed: 0.0 _TIME_


> `(bowtie2 -p 16 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_PRO-seq -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/prealignments/K562_PRO-seq_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
484143694 reads; of these:
  484143694 (100.00%) were unpaired; of these:
    439463330 (90.77%) aligned 0 times
    44680364 (9.23%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.23% overall alignment rate

> `Aligned_reads_human_rDNA`	44680364.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.23	PEPPRO	_RES_

### Map to genome (02-27 11:49:16) elapsed: 2987.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam`  

> `bowtie2 -p 16 --very-sensitive -X 2000 --rg-id K562_PRO-seq -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/prealignments/K562_PRO-seq_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/tmpn34_mq1l -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam` (306773,306779,306780)
<pre>
439463330 reads; of these:
  439463330 (100.00%) were unpaired; of these:
    5393540 (1.23%) aligned 0 times
    328159166 (74.67%) aligned exactly 1 time
    105910624 (24.10%) aligned >1 times
98.77% overall alignment rate
[bam_sort_core] merging from 131 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 4:01:46. Running peak memory: 3.842GB.  
  PID: 306779;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 306773;	Command: bowtie2;	Return code: 0;	Memory used: 3.842GB  
  PID: 306780;	Command: samtools;	Return code: 0;	Memory used: 0.883GB


> `samtools view -q 10 -b -@ 16 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_fail_qc.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam` (332405)
<pre>
</pre>
Command completed. Elapsed time: 0:13:37. Running peak memory: 3.842GB.  
  PID: 332405;	Command: samtools;	Return code: 0;	Memory used: 0.022GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	434069790	PEPPRO	_RES_

> `QC_filtered_reads`	46741044	PEPPRO	_RES_

> `Aligned_reads`	387328746	PEPPRO	_RES_

> `Alignment_rate`	80.0	PEPPRO	_RES_

> `Total_efficiency`	78.0	PEPPRO	_RES_

> `Read_depth`	29.61	PEPPRO	_RES_

### Compress all unmapped read files (02-27 17:01:58) elapsed: 18761.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/prealignments/K562_PRO-seq_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/prealignments/K562_PRO-seq_human_rDNA_unmap.fq` (339840)
<pre>
</pre>
Command completed. Elapsed time: 0:12:39. Running peak memory: 3.842GB.  
  PID: 339840;	Command: pigz;	Return code: 0;	Memory used: 0.012GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam` (341105)
<pre>
</pre>
Command completed. Elapsed time: 0:06:37. Running peak memory: 3.842GB.  
  PID: 341105;	Command: samtools;	Return code: 0;	Memory used: 0.027GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	9149071	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam` (342278)
<pre>
</pre>
Command completed. Elapsed time: 0:05:56. Running peak memory: 3.842GB.  
  PID: 342278;	Command: samtools;	Return code: 0;	Memory used: 0.027GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_noMT.bam` (342840,342845,342846,342847)
<pre>
</pre>
Command completed. Elapsed time: 0:07:44. Running peak memory: 3.842GB.  
  PID: 342845;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 342840;	Command: samtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 342846;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 342847;	Command: xargs;	Return code: 0;	Memory used: 0.161GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_noMT.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam` (343613)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.842GB.  
  PID: 343613;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam` (343623)
<pre>
</pre>
Command completed. Elapsed time: 0:05:55. Running peak memory: 3.842GB.  
  PID: 343623;	Command: samtools;	Return code: 0;	Memory used: 0.027GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-27 17:55:07) elapsed: 3190.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -c 16 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv` (345901)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/tmp_K562_PRO-seq_sort_htu017vb'
Processing with 16 cores...
Discarding 81 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1']
Keeping 114 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270411v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:07:00. Running peak memory: 12.017GB.  
  PID: 345901;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 12.017GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv`

> `NRF`	0.28	PEPPRO	_RES_

> `PBC1`	0.54	PEPPRO	_RES_

> `PBC2`	3.89	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_unmap.bam`  

> `samtools view -b -@ 16 -f 4  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_unmap.bam` (346626)
<pre>
</pre>
Command completed. Elapsed time: 0:01:17. Running peak memory: 12.017GB.  
  PID: 346626;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools view -c -f 4 -@ 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam`

> `Unmapped_reads`	5393540	PEPPRO	_RES_

### Split BAM by strand (02-27 18:04:36) elapsed: 569.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam` (346855)
<pre>
</pre>
Command completed. Elapsed time: 0:30:44. Running peak memory: 12.017GB.  
  PID: 346855;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam` (350436)
<pre>
</pre>
Command completed. Elapsed time: 0:31:47. Running peak memory: 12.017GB.  
  PID: 350436;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-27 19:07:06) elapsed: 3751.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (353879)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.017GB.  
  PID: 353879;	Command: sed;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/plus_TSS.tsv -p ends -c 16 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_plus_TssEnrichment.txt` (353885)
<pre>
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 12.017GB.  
  PID: 353885;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.268GB


> `TSS_coding_score`	13.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/minus_TSS.tsv -p ends -c 16 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_minus_TssEnrichment.txt` (353972)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 12.017GB.  
  PID: 353972;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.778GB


> `TSS_non-coding_score`	5.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_minus_TssEnrichment.txt` (354049)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 12.017GB.  
  PID: 354049;	Command: Rscript;	Return code: 0;	Memory used: 0.072GB

> `TSS enrichment`	QC_hg38/K562_PRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt` (354077,354078,354079,354080)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.017GB.  
  PID: 354077;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 354079;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 354078;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 354080;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt` (354082)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.017GB.  
  PID: 354082;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-27 19:08:45) elapsed: 98.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_tss.bed` (354084,354085)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 12.017GB.  
  PID: 354084;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 354085;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed` (354089,354090)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.017GB.  
  PID: 354089;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 354090;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_TSS_density.bed` (354092,354093,354094,354095)
<pre>
</pre>
Command completed. Elapsed time: 0:10:19. Running peak memory: 12.017GB.  
  PID: 354093;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 354095;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 354092;	Command: bedtools;	Return code: 0;	Memory used: 0.103GB  
  PID: 354094;	Command: sort;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_gene_body_density.bed` (355145,355150,355151)
<pre>
</pre>
Command completed. Elapsed time: 0:17:54. Running peak memory: 12.017GB.  
  PID: 355145;	Command: bedtools;	Return code: 0;	Memory used: 1.019GB  
  PID: 355151;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 355150;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed` (357194,357201,357202)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.017GB.  
  PID: 357194;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 357202;	Command: env;	Return code: 0;	Memory used: 0.005GB  
  PID: 357201;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed` (357207)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 12.017GB.  
  PID: 357207;	Command: Rscript;	Return code: 0;	Memory used: 0.215GB

> `Pause index`	QC_hg38/K562_PRO-seq_pause_index.pdf	Pause index	QC_hg38/K562_PRO-seq_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed.gz`  

> `pigz -f -p 16 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed` (357233)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.017GB.  
  PID: 357233;	Command: pigz;	Return code: 0;	Memory used: 0.004GB


### Calculate Fraction of Reads in pre-mature mRNA (02-27 19:37:07) elapsed: 1702.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam`
387328746 138157932

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam`
387328746 131399421

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/hg38_gene_sort.bed` (358495,358502)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 12.017GB.  
  PID: 358495;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 358502;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_gene_coverage.bed` (358504)
<pre>
</pre>
Command completed. Elapsed time: 0:17:47. Running peak memory: 12.017GB.  
  PID: 358504;	Command: bedtools;	Return code: 0;	Memory used: 1.018GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed.gz` (360479)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.017GB.  
  PID: 360479;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 16 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed` (360486)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.017GB.  
  PID: 360486;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (02-27 20:06:22) elapsed: 1755.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed` (360494)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 12.017GB.  
  PID: 360494;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR"` (360497)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.017GB.  
  PID: 360497;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR_sort.bed` (360498,360499,360500,360501)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.017GB.  
  PID: 360498;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 360499;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 360501;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 360500;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_plus_coverage.bed` (360504)
<pre>
</pre>
Command completed. Elapsed time: 0:05:40. Running peak memory: 12.017GB.  
  PID: 360504;	Command: bedtools;	Return code: 0;	Memory used: 0.104GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_minus_coverage.bed` (361141)
<pre>
</pre>
Command completed. Elapsed time: 0:04:50. Running peak memory: 12.017GB.  
  PID: 361141;	Command: bedtools;	Return code: 0;	Memory used: 0.163GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR"` (361759)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.017GB.  
  PID: 361759;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR_sort.bed` (361760,361761,361762,361763)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.017GB.  
  PID: 361760;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 361761;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 361763;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 361762;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_plus_coverage.bed` (361765)
<pre>
</pre>
Command completed. Elapsed time: 0:04:30. Running peak memory: 12.017GB.  
  PID: 361765;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_minus_coverage.bed` (362262)
<pre>
</pre>
Command completed. Elapsed time: 0:04:58. Running peak memory: 12.017GB.  
  PID: 362262;	Command: bedtools;	Return code: 0;	Memory used: 0.064GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer_sort.bed` (362763,362764,362765,362766)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.017GB.  
  PID: 362763;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 362764;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 362766;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 362765;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_plus_coverage.bed` (362769)
<pre>
</pre>
Command completed. Elapsed time: 0:05:02. Running peak memory: 12.017GB.  
  PID: 362769;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_minus_coverage.bed` (363315)
<pre>
</pre>
Command completed. Elapsed time: 0:04:20. Running peak memory: 12.017GB.  
  PID: 363315;	Command: bedtools;	Return code: 0;	Memory used: 0.068GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Exon_sort.bed` (363812,363813,363814,363815)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 12.017GB.  
  PID: 363812;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 363813;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 363815;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 363814;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_plus_coverage.bed` (363820)
<pre>
</pre>
Command completed. Elapsed time: 0:05:58. Running peak memory: 12.017GB.  
  PID: 363820;	Command: bedtools;	Return code: 0;	Memory used: 0.466GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_minus_coverage.bed` (364404)
<pre>
</pre>
Command completed. Elapsed time: 0:05:52. Running peak memory: 12.017GB.  
  PID: 364404;	Command: bedtools;	Return code: 0;	Memory used: 0.186GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Intron_sort.bed` (365025,365031,365032,365033)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 12.017GB.  
  PID: 365025;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 365032;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 365031;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 365033;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_plus_coverage.bed` (365037)
<pre>
</pre>
Command completed. Elapsed time: 0:07:50. Running peak memory: 12.017GB.  
  PID: 365037;	Command: bedtools;	Return code: 0;	Memory used: 0.45GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_minus_coverage.bed` (366060)
<pre>
</pre>
Command completed. Elapsed time: 0:07:58. Running peak memory: 12.017GB.  
  PID: 366060;	Command: bedtools;	Return code: 0;	Memory used: 0.433GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_sort.bed` (366806,366812,366813,366814)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.017GB.  
  PID: 366806;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 366812;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 366814;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 366813;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_plus_coverage.bed` (366817)
<pre>
</pre>
Command completed. Elapsed time: 0:05:41. Running peak memory: 12.017GB.  
  PID: 366817;	Command: bedtools;	Return code: 0;	Memory used: 0.581GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_minus_coverage.bed` (367390)
<pre>
</pre>
Command completed. Elapsed time: 0:05:36. Running peak memory: 12.017GB.  
  PID: 367390;	Command: bedtools;	Return code: 0;	Memory used: 0.471GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region"` (367963)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.017GB.  
  PID: 367963;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed` (367968,367969,367970,367971)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.017GB.  
  PID: 367968;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 367970;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 367969;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 367971;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (367974)
<pre>
</pre>
Command completed. Elapsed time: 0:05:05. Running peak memory: 12.017GB.  
  PID: 367974;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_minus_coverage.bed` (368543)
<pre>
</pre>
Command completed. Elapsed time: 0:05:25. Running peak memory: 12.017GB.  
  PID: 368543;	Command: bedtools;	Return code: 0;	Memory used: 0.165GB


### Plot cFRiF/FRiF (02-27 21:25:19) elapsed: 4737.0 _TIME_


> `samtools view -@ 16 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_PRO-seq -z 3099922541 -n 192750289 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (369468)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:48. Running peak memory: 12.017GB.  
  PID: 369468;	Command: Rscript;	Return code: 0;	Memory used: 0.803GB

> `cFRiF`	QC_hg38/K562_PRO-seq_cFRiF.pdf	cFRiF	QC_hg38/K562_PRO-seq_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_PRO-seq -z 3099922541 -n 192750289 -y frif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (369529)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 12.017GB.  
  PID: 369529;	Command: Rscript;	Return code: 0;	Memory used: 0.502GB

> `FRiF`	QC_hg38/K562_PRO-seq_FRiF.pdf	FRiF	QC_hg38/K562_PRO-seq_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-27 21:27:09) elapsed: 110.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/hg38_exons_sort.bed` (369568,369569)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 12.017GB.  
  PID: 369569;	Command: bedtools;	Return code: 0;	Memory used: 0.08GB  
  PID: 369568;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/hg38_introns_sort.bed` (369576,369577,369578)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 12.017GB.  
  PID: 369576;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 369578;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 369577;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_coverage.bed` (369584)
<pre>
</pre>
Command completed. Elapsed time: 0:10:43. Running peak memory: 12.017GB.  
  PID: 369584;	Command: bedtools;	Return code: 0;	Memory used: 0.291GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_coverage.bed` (370713)
<pre>
</pre>
Command completed. Elapsed time: 0:17:26. Running peak memory: 12.017GB.  
  PID: 370713;	Command: bedtools;	Return code: 0;	Memory used: 0.486GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/387.328746)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_rpkm.bed` (372653,372659,372660)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.017GB.  
  PID: 372653;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 372660;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 372659;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/387.328746)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_rpkm.bed` (372662,372663,372664)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.017GB.  
  PID: 372662;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 372664;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 372663;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_rpkm.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed` (372667,372669,372670)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.017GB.  
  PID: 372667;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 372670;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 372669;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed --annotate` (372676)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 12.017GB.  
  PID: 372676;	Command: Rscript;	Return code: 0;	Memory used: 0.27GB

> `mRNA contamination`	QC_hg38/K562_PRO-seq_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_PRO-seq_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed.gz`  

> `pigz -f -p 16 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed` (372731)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.017GB.  
  PID: 372731;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (02-27 21:55:37) elapsed: 1707.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam` (372740)
<pre>
</pre>
Command completed. Elapsed time: 0:03:01. Running peak memory: 12.017GB.  
  PID: 372740;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_smooth_body_0-mer.bw -p 10 --variable-step --tail-edge` (372925)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_plus_cuttrace_rxd14k9c'
Processing with 5 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 90 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1']
Keeping 105 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 105 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_exact_body_0-mer.bw'
Merging 105 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:13:54. Running peak memory: 12.017GB.  
  PID: 372925;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.595GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam` (375569)
<pre>
</pre>
Command completed. Elapsed time: 0:02:57. Running peak memory: 12.017GB.  
  PID: 375569;	Command: samtools;	Return code: 0;	Memory used: 0.02GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_smooth_body_0-mer.bw -p 10 --variable-step --tail-edge` (375945)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_minus_cuttrace_pjoixdgn'
Processing with 5 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 93 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 102 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270411v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270336v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 102 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_exact_body_0-mer.bw'
Merging 102 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:13:02. Running peak memory: 12.017GB.  
  PID: 375945;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.835GB

Starting cleanup: 52 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  13:05:37
*  Total elapsed time (all runs):  17:52:40
*         Peak memory (this run):  12.0173 GB
*        Pipeline completed time: 2020-02-27 22:28:40
