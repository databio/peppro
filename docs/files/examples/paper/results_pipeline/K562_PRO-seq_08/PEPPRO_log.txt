### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_08 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_8pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty`
*         Compute host:  udc-aw29-26a
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/
*  Pipeline started at:   (06-11 17:02:38) elapsed: 1.0 _TIME_

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
*              `cores`:  `12`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_8pct.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `12000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `False`
*        `sample_name`:  `K562_PRO-seq_08`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_8pct.fastq.gz

> `File_mb`	2954.27	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 17:02:38) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/K562_PRO-seq_08.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_PRO_8pct.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/K562_PRO-seq_08.fastq.gz` (368927)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 368927;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/K562_PRO-seq_08.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/K562_PRO-seq_08_R1.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/K562_PRO-seq_08.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/K562_PRO-seq_08_R1.fastq` (368928)
<pre>
</pre>
Command completed. Elapsed time: 0:02:14. Running peak memory: 0.002GB.  
  PID: 368928;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	39716621	PEPPRO	_RES_

> `Fastq_reads`	39716621	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/K562_PRO-seq_08.fastq.gz']

### FASTQ processing:  (06-11 17:06:01) elapsed: 203.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/K562_PRO-seq_08_R1_processed.fastq`  

> `(cutadapt -j 12 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/K562_PRO-seq_08_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/K562_PRO-seq_08_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/K562_PRO-seq_08_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt/K562_PRO-seq_08_R1_cutadapt.txt` (369343)
<pre>
</pre>
Command completed. Elapsed time: 0:02:00. Running peak memory: 4.208GB.  
  PID: 369343;	Command: cutadapt;	Return code: 0;	Memory used: 4.208GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/K562_PRO-seq_08_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/K562_PRO-seq_08_R1_processed.fastq` (369493,369494)
<pre>
</pre>
Command completed. Elapsed time: 0:01:50. Running peak memory: 4.208GB.  
  PID: 369493;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 369494;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt/K562_PRO-seq_08_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	33831973.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt/K562_PRO-seq_08_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt/K562_PRO-seq_08_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/K562_PRO-seq_08_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	995532.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	2.5066	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastqc/K562_PRO-seq_08_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (369975)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.208GB.  
  PID: 369975;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	38721089	PEPPRO	_RES_

> `Trim_loss_rate`	2.51	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/K562_PRO-seq_08_R1_processed.fastq` (369998)
<pre>
Started analysis of K562_PRO-seq_08_R1_processed.fastq
Approx 5% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 10% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 15% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 20% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 25% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 30% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 35% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 40% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 45% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 50% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 55% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 60% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 65% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 70% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 75% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 80% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 85% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 90% complete for K562_PRO-seq_08_R1_processed.fastq
Approx 95% complete for K562_PRO-seq_08_R1_processed.fastq
Analysis complete for K562_PRO-seq_08_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:02:29. Running peak memory: 4.208GB.  
  PID: 369998;	Command: fastqc;	Return code: 0;	Memory used: 0.204GB

> `FastQC report r1`	fastqc/K562_PRO-seq_08_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/processed_R1.flag` (370147)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.208GB.  
  PID: 370147;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Plot adapter insertion distribution (06-11 17:13:15) elapsed: 433.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt/K562_PRO-seq_08_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt/K562_PRO-seq_08_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt` (370148)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 4.208GB.  
  PID: 370148;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_08_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_08_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt/K562_PRO-seq_08_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 17:13:21) elapsed: 7.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt/K562_PRO-seq_08_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt/K562_PRO-seq_08_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt/K562_PRO-seq_08_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt/K562_PRO-seq_08_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt/K562_PRO-seq_08_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2313	PEPPRO	_RES_

### Prealignments (06-11 17:13:21) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 17:13:21) elapsed: 0.0 _TIME_


> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_PRO-seq_08 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/K562_PRO-seq_08_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/prealignments/K562_PRO-seq_08_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
38721089 reads; of these:
  38721089 (100.00%) were unpaired; of these:
    35147182 (90.77%) aligned 0 times
    3573907 (9.23%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.23% overall alignment rate

> `Aligned_reads_human_rDNA`	3573907.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.23	PEPPRO	_RES_

### Map to genome (06-11 17:18:02) elapsed: 280.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam`  

> `bowtie2 -p 12 --very-sensitive --rg-id K562_PRO-seq_08 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/prealignments/K562_PRO-seq_08_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/tmp4z0_s6qa -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_temp.bam` (370659,370660,370661)
<pre>
35147182 reads; of these:
  35147182 (100.00%) were unpaired; of these:
    432298 (1.23%) aligned 0 times
    26245142 (74.67%) aligned exactly 1 time
    8469742 (24.10%) aligned >1 times
98.77% overall alignment rate
[bam_sort_core] merging from 11 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:27:31. Running peak memory: 4.208GB.  
  PID: 370660;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 370659;	Command: bowtie2;	Return code: 0;	Memory used: 3.717GB  
  PID: 370661;	Command: samtools;	Return code: 0;	Memory used: 0.889GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam` (374389)
<pre>
</pre>
Command completed. Elapsed time: 0:01:08. Running peak memory: 4.208GB.  
  PID: 374389;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	34714884	PEPPRO	_RES_

> `QC_filtered_reads`	3735729	PEPPRO	_RES_

> `Aligned_reads`	30979155	PEPPRO	_RES_

> `Alignment_rate`	80.01	PEPPRO	_RES_

> `Total_efficiency`	78.0	PEPPRO	_RES_

> `Read_depth`	3.98	PEPPRO	_RES_

### Compress all unmapped read files (06-11 18:04:13) elapsed: 2772.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/prealignments/K562_PRO-seq_08_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/prealignments/K562_PRO-seq_08_human_rDNA_unmap.fq` (377456)
<pre>
</pre>
Command completed. Elapsed time: 0:01:28. Running peak memory: 4.208GB.  
  PID: 377456;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_temp.bam` (377848)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 4.208GB.  
  PID: 377848;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	731070	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam` (377921)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 4.208GB.  
  PID: 377921;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/chr_sizes.bed` (377979,377980,377981,377982)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.208GB.  
  PID: 377979;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 377981;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 377980;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 377982;	Command: grep;	Return code: 0;	Memory used: 0.002GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_noMT.bam` (377985)
<pre>
</pre>
Command completed. Elapsed time: 0:00:44. Running peak memory: 4.208GB.  
  PID: 377985;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam` (378096)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.208GB.  
  PID: 378096;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam` (378097)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 4.208GB.  
  PID: 378097;	Command: samtools;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-11 18:09:49) elapsed: 335.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_bamQC.tsv` (378448)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/tmp_K562_PRO-seq_08_sort_fdwp9579'
Processing with 12 cores...
Discarding 100 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 95 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 4.208GB.  
  PID: 378448;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.733GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_bamQC.tsv`

> `NRF`	0.84	PEPPRO	_RES_

> `PBC1`	0.92	PEPPRO	_RES_

> `PBC2`	16.03	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_unmap.bam`  

> `samtools view -b -@ 12 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_unmap.bam` (378897)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 4.208GB.  
  PID: 378897;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_temp.bam`

> `Unmapped_reads`	432298	PEPPRO	_RES_

### Split BAM by strand (06-11 18:10:40) elapsed: 52.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam` (378950)
<pre>
</pre>
Command completed. Elapsed time: 0:02:38. Running peak memory: 4.208GB.  
  PID: 378950;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam` (379393)
<pre>
</pre>
Command completed. Elapsed time: 0:02:32. Running peak memory: 4.208GB.  
  PID: 379393;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-11 18:15:51) elapsed: 310.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (380113)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.208GB.  
  PID: 380113;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_plus_TssEnrichment.txt` (380114)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 4.208GB.  
  PID: 380114;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.517GB


> `TSS_coding_score`	13.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_minus_TssEnrichment.txt` (380155)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 4.208GB.  
  PID: 380155;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.664GB


> `TSS_non-coding_score`	5.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_minus_TssEnrichment.txt` (380193)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 4.208GB.  
  PID: 380193;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `TSS enrichment`	QC_hg38/K562_PRO-seq_08_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_08_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt` (380227,380228,380229,380230)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.208GB.  
  PID: 380227;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 380229;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 380228;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 380230;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_keep.txt` (380232)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.208GB.  
  PID: 380232;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-11 18:16:13) elapsed: 22.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_ensembl_tss.bed` (380234,380235)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 4.208GB.  
  PID: 380234;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 380235;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_ensembl_gene_body.bed` (380241,380242)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.208GB.  
  PID: 380241;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 380242;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_TSS_density.bed` (380244,380245,380246,380247)
<pre>
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 4.208GB.  
  PID: 380247;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 380244;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB  
  PID: 380246;	Command: sort;	Return code: 0;	Memory used: 0.008GB  
  PID: 380245;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_gene_body_density.bed` (380341,380342,380343)
<pre>
</pre>
Command completed. Elapsed time: 0:01:06. Running peak memory: 4.208GB.  
  PID: 380342;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 380341;	Command: bedtools;	Return code: 0;	Memory used: 0.1GB  
  PID: 380343;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_pause_index.bed`  
awk: cmd. line:1: { if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}
awk: cmd. line:1:             ^ syntax error
awk: cmd. line:1: { if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}
awk: cmd. line:1:                                                                                                                      ^ syntax error

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/tmpc39o6tyo` (380613,380614,380615)
<pre>
</pre>
Got SIGTERM. Failing gracefully... (06-11 18:55:13) elapsed: 2340.0 _TIME_
Child process 380613 (join) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/recover.lock.QC_hg38__K562_PRO-seq_08_pause_index.bed

### Pipeline failed at:  (06-11 18:55:13) elapsed: 0.0 _TIME_

Total time: 1:52:36
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_08 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_8pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba26-34c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/
*  Pipeline started at:   (06-11 19:03:07) elapsed: 1.0 _TIME_

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
*              `cores`:  `12`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_8pct.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `12000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `K562_PRO-seq_08`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_8pct.fastq.gz

> `File_mb`	2954.27	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 19:03:07) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/K562_PRO-seq_08.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/K562_PRO-seq_08.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/K562_PRO-seq_08_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/K562_PRO-seq_08.fastq.gz']

### FASTQ processing:  (06-11 19:03:07) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-11 19:03:07) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt/K562_PRO-seq_08_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_08_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_08_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-11 19:03:07) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 19:03:07) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/prealignments/K562_PRO-seq_08_human_rDNA_unmap.fq

### Map to genome (06-11 19:03:07) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam`  

### Compress all unmapped read files (06-11 19:03:07) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/prealignments/K562_PRO-seq_08_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-11 19:03:07) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_unmap.bam`  

### Split BAM by strand (06-11 19:03:07) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam`  

### Calculate TSS enrichment (06-11 19:03:07) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/K562_PRO-seq_08_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_08_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-11 19:03:07) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_gene_body_density.bed`  
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/lock.QC_hg38__K562_PRO-seq_08_pause_index.bed
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/tmptrrv1kh3` (198993,198994,198995)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 198993;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 198995;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 198994;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/tmptrrv1kh3 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-11 19:03:08) elapsed: 0.0 _TIME_

Total time: 0:00:03
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/tmptrrv1kh3 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_08 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_8pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aw29-25b
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/
*  Pipeline started at:   (06-14 21:09:09) elapsed: 5.0 _TIME_

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
*              `cores`:  `12`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_8pct.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `12000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `K562_PRO-seq_08`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_8pct.fastq.gz

> `File_mb`	2954.27	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-14 21:09:09) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/K562_PRO-seq_08.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/K562_PRO-seq_08.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/K562_PRO-seq_08_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/K562_PRO-seq_08.fastq.gz']

### FASTQ processing:  (06-14 21:09:09) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-14 21:09:09) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/cutadapt/K562_PRO-seq_08_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_08_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_08_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-14 21:09:09) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:09:09) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/prealignments/K562_PRO-seq_08_human_rDNA_unmap.fq

### Map to genome (06-14 21:09:09) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam`  

### Compress all unmapped read files (06-14 21:09:09) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/prealignments/K562_PRO-seq_08_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:09:09) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_unmap.bam`  

### Split BAM by strand (06-14 21:09:09) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam`  

### Calculate TSS enrichment (06-14 21:09:09) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/K562_PRO-seq_08_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_08_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:09:09) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/tmplwfed8j9` (378607,378610,378612)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.006GB.  
  PID: 378607;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 378612;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 378610;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/tmplwfed8j9 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0267261) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/tmplwfed8j9 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_pause_index.bed` (378641)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.006GB.  
  PID: 378641;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.59	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_pause_index.bed` (378649)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.285GB.  
  PID: 378649;	Command: Rscript;	Return code: 0;	Memory used: 0.285GB

> `Pause index`	QC_hg38/K562_PRO-seq_08_pause_index.pdf	Pause index	QC_hg38/K562_PRO-seq_08_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_pause_index.bed` (378688)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.285GB.  
  PID: 378688;	Command: pigz;	Return code: 0;	Memory used: 0.004GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:09:15) elapsed: 6.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam`
30979155 11051521

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam`
30979155 10511222

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/signal_hg38/K562_PRO-seq_08_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_gene_sort.bed` (379102,379103)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.285GB.  
  PID: 379102;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 379103;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/signal_hg38/K562_PRO-seq_08_gene_coverage.bed` (379106)
<pre>
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 0.285GB.  
  PID: 379106;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/hg38_annotations.bed.gz` (379171)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.285GB.  
  PID: 379171;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/hg38_annotations.bed` (379174)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.285GB.  
  PID: 379174;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:11:14) elapsed: 118.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/raw/hg38_annotations.bed` (379183)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.285GB.  
  PID: 379183;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Enhancer_sort.bed` (379187,379188,379189,379190)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.285GB.  
  PID: 379187;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 379188;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 379190;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 379189;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Enhancer_plus_coverage.bed` (379192)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 0.285GB.  
  PID: 379192;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Enhancer_minus_coverage.bed` (379223)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 0.285GB.  
  PID: 379223;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter_sort.bed` (379368,379369,379370,379371)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.285GB.  
  PID: 379368;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 379369;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 379371;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 379370;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Promoter_plus_coverage.bed` (379373)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 0.285GB.  
  PID: 379373;	Command: bedtools;	Return code: 0;	Memory used: 0.027GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Promoter_minus_coverage.bed` (379397)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 0.285GB.  
  PID: 379397;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter_Flanking_Region"` (379435)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.285GB.  
  PID: 379435;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter_Flanking_Region_sort.bed` (379436,379437,379438,379439)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.285GB.  
  PID: 379436;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 379438;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 379437;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 379439;	Command: bedtools;	Return code: 0;	Memory used: 0.05GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Promoter_Flanking_Region_plus_coverage.bed` (379442)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 0.285GB.  
  PID: 379442;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Promoter_Flanking_Region_minus_coverage.bed` (386555)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 0.285GB.  
  PID: 386555;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/5_UTR"` (390479)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.285GB.  
  PID: 390479;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/5_UTR_sort.bed` (390482,390484,390486,390488)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.285GB.  
  PID: 390482;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 390484;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 390488;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 390486;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_5_UTR_plus_coverage.bed` (390498)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 0.285GB.  
  PID: 390498;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_5_UTR_minus_coverage.bed` (397012)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 0.285GB.  
  PID: 397012;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/3_UTR"` (407030)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.285GB.  
  PID: 407030;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/3_UTR_sort.bed` (407031,407032,407033,407034)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.285GB.  
  PID: 407031;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 407032;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 407034;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB  
  PID: 407033;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_3_UTR_plus_coverage.bed` (407036)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 0.285GB.  
  PID: 407036;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_3_UTR_minus_coverage.bed` (407088)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 0.285GB.  
  PID: 407088;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Exon_sort.bed` (407319,407320,407321,407322)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.285GB.  
  PID: 407319;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 407320;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 407322;	Command: bedtools;	Return code: 0;	Memory used: 0.17GB  
  PID: 407321;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Exon_plus_coverage.bed` (407327)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 0.285GB.  
  PID: 407327;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Exon_minus_coverage.bed` (407365)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 0.285GB.  
  PID: 407365;	Command: bedtools;	Return code: 0;	Memory used: 0.023GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Intron_sort.bed` (407385,407386,407387,407388)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.285GB.  
  PID: 407385;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 407387;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 407386;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 407388;	Command: bedtools;	Return code: 0;	Memory used: 0.083GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Intron_plus_coverage.bed` (407393)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 0.285GB.  
  PID: 407393;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Intron_minus_coverage.bed` (407440)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 0.285GB.  
  PID: 407440;	Command: bedtools;	Return code: 0;	Memory used: 0.047GB


### Plot cFRiF/FRiF (06-14 21:16:49) elapsed: 335.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_08 -z 3099922541 -n 15416586 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Intron_plus_coverage.bed` (407484)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 0.465GB.  
  PID: 407484;	Command: Rscript;	Return code: 0;	Memory used: 0.465GB

> `cFRiF`	QC_hg38/K562_PRO-seq_08_cFRiF.pdf	cFRiF	QC_hg38/K562_PRO-seq_08_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_08 -z 3099922541 -n 15416586 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_Intron_plus_coverage.bed` (413036)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 0.51GB.  
  PID: 413036;	Command: Rscript;	Return code: 0;	Memory used: 0.51GB

> `FRiF`	QC_hg38/K562_PRO-seq_08_FRiF.pdf	FRiF	QC_hg38/K562_PRO-seq_08_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 21:17:52) elapsed: 63.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_exons_sort.bed` (413071,413072)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.51GB.  
  PID: 413072;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB  
  PID: 413071;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_introns_sort.bed` (413078,413079,413080)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.51GB.  
  PID: 413078;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 413080;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 413079;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_exons_coverage.bed` (413088)
<pre>
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 0.51GB.  
  PID: 413088;	Command: bedtools;	Return code: 0;	Memory used: 0.031GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_introns_coverage.bed` (429644)
<pre>
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 0.51GB.  
  PID: 429644;	Command: bedtools;	Return code: 0;	Memory used: 0.053GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/30.979155)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_exons_rpkm.bed` (429726,429727,429728)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.51GB.  
  PID: 429726;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 429728;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 429727;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/30.979155)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_introns_rpkm.bed` (429730,429731,429732)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.51GB.  
  PID: 429730;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 429732;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 429731;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_exon_intron_ratios.bed` (429735,429736,429737)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.51GB.  
  PID: 429735;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 429737;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 429736;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.41	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_exon_intron_ratios.bed --annotate` (429743)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.51GB.  
  PID: 429743;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/K562_PRO-seq_08_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_PRO-seq_08_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/QC_hg38/K562_PRO-seq_08_exon_intron_ratios.bed` (429768)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.51GB.  
  PID: 429768;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (06-14 21:19:43) elapsed: 111.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/signal_hg38/K562_PRO-seq_08_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/signal_hg38/K562_PRO-seq_08_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam` (429776)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 0.51GB.  
  PID: 429776;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/signal_hg38/K562_PRO-seq_08_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/signal_hg38/K562_PRO-seq_08_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 30979155.0` (429790)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_plus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_08_plus_cuttrace_yx9e42iv'
Processing with 4 cores...
Discarding 113 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr22_KI270732v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrEBV']
Keeping 82 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 82 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/signal_hg38/K562_PRO-seq_08_plus_exact_body_0-mer.bw'
Merging 82 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/signal_hg38/K562_PRO-seq_08_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:37. Running peak memory: 2.631GB.  
  PID: 429790;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.631GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/signal_hg38/K562_PRO-seq_08_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/signal_hg38/K562_PRO-seq_08_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam` (18142)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 2.631GB.  
  PID: 18142;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/signal_hg38/K562_PRO-seq_08_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/signal_hg38/K562_PRO-seq_08_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 30979155.0` (18179)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/aligned_hg38/K562_PRO-seq_08_minus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_08_minus_cuttrace_7tl2zsla'
Processing with 4 cores...
stdin is empty of data
Discarding 112 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1']
Keeping 83 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/signal_hg38/K562_PRO-seq_08_minus_exact_body_0-mer.bw'
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_08/signal_hg38/K562_PRO-seq_08_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:34. Running peak memory: 2.631GB.  
  PID: 18179;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.548GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:26:19
*  Total elapsed time (all runs):  1:59:46
*         Peak memory (this run):  2.6308 GB
*        Pipeline completed time: 2020-06-14 21:35:23
