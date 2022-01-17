### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_60 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_60pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty`
*         Compute host:  udc-ba26-16
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/
*  Pipeline started at:   (06-11 17:15:22) elapsed: 6.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_60pctRNA.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_60`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_60pctRNA.fastq.gz

> `File_mb`	5477.44	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 17:15:23) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/K562_RNA-seq_60.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_60pctRNA.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/K562_RNA-seq_60.fastq.gz` (47287)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 47287;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/K562_RNA-seq_60.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/K562_RNA-seq_60_R1.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/K562_RNA-seq_60.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/K562_RNA-seq_60_R1.fastq` (47288)
<pre>
</pre>
Command completed. Elapsed time: 0:02:46. Running peak memory: 0.002GB.  
  PID: 47288;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	70000000	PEPPRO	_RES_

> `Fastq_reads`	70000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/K562_RNA-seq_60.fastq.gz']

### FASTQ processing:  (06-11 17:19:51) elapsed: 268.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/K562_RNA-seq_60_R1_processed.fastq`  

> `(cutadapt -j 12 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/K562_RNA-seq_60_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/K562_RNA-seq_60_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/K562_RNA-seq_60_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt/K562_RNA-seq_60_R1_cutadapt.txt` (48353)
<pre>
</pre>
Command completed. Elapsed time: 0:02:27. Running peak memory: 3.995GB.  
  PID: 48353;	Command: cutadapt;	Return code: 0;	Memory used: 3.995GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/K562_RNA-seq_60_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/K562_RNA-seq_60_R1_processed.fastq` (49360,49361)
<pre>
</pre>
Command completed. Elapsed time: 0:02:08. Running peak memory: 3.995GB.  
  PID: 49361;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB  
  PID: 49360;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt/K562_RNA-seq_60_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	38768961.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt/K562_RNA-seq_60_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt/K562_RNA-seq_60_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/K562_RNA-seq_60_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	700276.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	1.0004	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastqc/K562_RNA-seq_60_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (50334)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.995GB.  
  PID: 50334;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	69299724	PEPPRO	_RES_

> `Trim_loss_rate`	1.0	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/K562_RNA-seq_60_R1_processed.fastq` (50455)
<pre>
Started analysis of K562_RNA-seq_60_R1_processed.fastq
Approx 5% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 10% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 15% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 20% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 25% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 30% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 35% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 40% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 45% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 50% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 55% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 60% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 65% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 70% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 75% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 80% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 85% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 90% complete for K562_RNA-seq_60_R1_processed.fastq
Approx 95% complete for K562_RNA-seq_60_R1_processed.fastq
Analysis complete for K562_RNA-seq_60_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:04:54. Running peak memory: 3.995GB.  
  PID: 50455;	Command: fastqc;	Return code: 0;	Memory used: 0.248GB

> `FastQC report r1`	fastqc/K562_RNA-seq_60_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/processed_R1.flag` (51064)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.995GB.  
  PID: 51064;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-11 17:30:39) elapsed: 648.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt/K562_RNA-seq_60_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt/K562_RNA-seq_60_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt` (51065)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.995GB.  
  PID: 51065;	Command: Rscript;	Return code: 0;	Memory used: 0.123GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_60_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_60_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt/K562_RNA-seq_60_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 17:30:48) elapsed: 9.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt/K562_RNA-seq_60_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt/K562_RNA-seq_60_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt/K562_RNA-seq_60_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt/K562_RNA-seq_60_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt/K562_RNA-seq_60_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2309	PEPPRO	_RES_

### Prealignments (06-11 17:30:48) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 17:30:48) elapsed: 0.0 _TIME_


> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_60 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/K562_RNA-seq_60_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/prealignments/K562_RNA-seq_60_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
69299724 reads; of these:
  69299724 (100.00%) were unpaired; of these:
    65807135 (94.96%) aligned 0 times
    3492589 (5.04%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
5.04% overall alignment rate

> `Aligned_reads_human_rDNA`	3492589.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	5.04	PEPPRO	_RES_

### Map to genome (06-11 17:39:53) elapsed: 545.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam`  

> `bowtie2 -p 12 --very-sensitive --rg-id K562_RNA-seq_60 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/prealignments/K562_RNA-seq_60_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/tmpm17w1gg4 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_temp.bam` (51822,51829,51830)
<pre>
65807135 reads; of these:
  65807135 (100.00%) were unpaired; of these:
    7219211 (10.97%) aligned 0 times
    38663688 (58.75%) aligned exactly 1 time
    19924236 (30.28%) aligned >1 times
89.03% overall alignment rate
[bam_sort_core] merging from 21 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:40:52. Running peak memory: 3.995GB.  
  PID: 51822;	Command: bowtie2;	Return code: 0;	Memory used: 3.716GB  
  PID: 51829;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 51830;	Command: samtools;	Return code: 0;	Memory used: 0.889GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam` (61025)
<pre>
</pre>
Command completed. Elapsed time: 0:03:13. Running peak memory: 3.995GB.  
  PID: 61025;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	58587924	PEPPRO	_RES_

> `QC_filtered_reads`	10244576	PEPPRO	_RES_

> `Aligned_reads`	48343348	PEPPRO	_RES_

> `Alignment_rate`	69.76	PEPPRO	_RES_

> `Total_efficiency`	69.06	PEPPRO	_RES_

> `Read_depth`	6.87	PEPPRO	_RES_

### Compress all unmapped read files (06-11 18:45:22) elapsed: 3929.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/prealignments/K562_RNA-seq_60_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/prealignments/K562_RNA-seq_60_human_rDNA_unmap.fq` (66626)
<pre>
</pre>
Command completed. Elapsed time: 0:02:42. Running peak memory: 3.995GB.  
  PID: 66626;	Command: pigz;	Return code: 0;	Memory used: 0.014GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_temp.bam` (67160)
<pre>
</pre>
Command completed. Elapsed time: 0:01:11. Running peak memory: 3.995GB.  
  PID: 67160;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	3031668	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam` (67369)
<pre>
</pre>
Command completed. Elapsed time: 0:00:52. Running peak memory: 3.995GB.  
  PID: 67369;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/chr_sizes.bed` (67838,67839,67840,67841)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.995GB.  
  PID: 67839;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 67841;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 67838;	Command: samtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 67840;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_noMT.bam` (67843)
<pre>
</pre>
Command completed. Elapsed time: 0:01:09. Running peak memory: 3.995GB.  
  PID: 67843;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam` (68057)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.995GB.  
  PID: 68057;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam` (68060)
<pre>
</pre>
Command completed. Elapsed time: 0:00:51. Running peak memory: 3.995GB.  
  PID: 68060;	Command: samtools;	Return code: 0;	Memory used: 0.015GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-11 18:54:49) elapsed: 567.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_bamQC.tsv` (68807)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/tmp_K562_RNA-seq_60_sort_n5un91ms'
Processing with 12 cores...
</pre>
Got SIGTERM. Failing gracefully... (06-11 18:55:15) elapsed: 26.0 _TIME_
Process ForkPoolWorker-12:
Traceback (most recent call last):
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/process.py", line 258, in _bootstrap
    self.run()
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/process.py", line 93, in run
    self._target(*self._args, **self._kwargs)
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/pool.py", line 108, in worker
    task = get()
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/queues.py", line 334, in get
    with self._rlock:
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/synchronize.py", line 96, in __enter__
    return self._semlock.__enter__()
KeyboardInterrupt
Process ForkPoolWorker-10:
Traceback (most recent call last):
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/process.py", line 258, in _bootstrap
    self.run()
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/process.py", line 93, in run
    self._target(*self._args, **self._kwargs)
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/pool.py", line 108, in worker
    task = get()
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/queues.py", line 334, in get
    with self._rlock:
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/synchronize.py", line 96, in __enter__
    return self._semlock.__enter__()
KeyboardInterrupt
Process ForkPoolWorker-8:
Traceback (most recent call last):
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/process.py", line 258, in _bootstrap
    self.run()
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/process.py", line 93, in run
    self._target(*self._args, **self._kwargs)
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/pool.py", line 108, in worker
    task = get()
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/queues.py", line 334, in get
    with self._rlock:
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/synchronize.py", line 96, in __enter__
    return self._semlock.__enter__()
KeyboardInterrupt
Process ForkPoolWorker-11:
Traceback (most recent call last):
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/process.py", line 258, in _bootstrap
    self.run()
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/process.py", line 93, in run
    self._target(*self._args, **self._kwargs)
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/pool.py", line 108, in worker
    task = get()
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/queues.py", line 335, in get
    res = self._reader.recv_bytes()
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/connection.py", line 216, in recv_bytes
    buf = self._recv_bytes(maxlength)
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/connection.py", line 407, in _recv_bytes
    buf = self._recv(4)
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/connection.py", line 379, in _recv
    chunk = read(handle, remaining)
KeyboardInterrupt
Process ForkPoolWorker-9:
Traceback (most recent call last):
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/process.py", line 258, in _bootstrap
    self.run()
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/process.py", line 93, in run
    self._target(*self._args, **self._kwargs)
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/pool.py", line 108, in worker
    task = get()
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/queues.py", line 334, in get
    with self._rlock:
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/synchronize.py", line 96, in __enter__
    return self._semlock.__enter__()
KeyboardInterrupt
Process ForkPoolWorker-5:
Traceback (most recent call last):
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/process.py", line 258, in _bootstrap
    self.run()
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/process.py", line 93, in run
    self._target(*self._args, **self._kwargs)
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/pool.py", line 108, in worker
    task = get()
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/queues.py", line 334, in get
    with self._rlock:
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/synchronize.py", line 96, in __enter__
    return self._semlock.__enter__()
KeyboardInterrupt
Process ForkPoolWorker-3:
Traceback (most recent call last):
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/process.py", line 258, in _bootstrap
    self.run()
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/process.py", line 93, in run
    self._target(*self._args, **self._kwargs)
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/pool.py", line 119, in worker
    result = (True, func(*args, **kwds))
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/pool.py", line 44, in mapstar
    return list(map(*args))
  File "/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py", line 148, in __call__
    readSE = getRead(self.fetch_chunk(chrom), isPE)
  File "/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py", line 94, in getRead
    readDict['query_name'].append(read.query_name)
KeyboardInterrupt
Traceback (most recent call last):
  File "/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py", line 249, in <module>
    good_chromosomes = qc.run()
  File "/home/jps3dp/.local/lib/python3.6/site-packages/pararead/processor.py", line 405, in run
    results = workers.map_async(self, nonempties).get(9999999)
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/pool.py", line 638, in get
    self.wait(timeout)
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/multiprocessing/pool.py", line 635, in wait
    self._event.wait(timeout)
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/threading.py", line 551, in wait
    signaled = self._cond.wait(timeout)
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/threading.py", line 299, in wait
    gotit = waiter.acquire(True, timeout)
KeyboardInterrupt
Error in atexit._run_exitfuncs:
Traceback (most recent call last):
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/shutil.py", line 480, in rmtree
    _rmtree_safe_fd(fd, path, onerror)
  File "/apps/software/standard/core/anaconda/5.2.0-py3.6/lib/python3.6/shutil.py", line 436, in _rmtree_safe_fd
    os.unlink(name, dir_fd=topfd)
KeyboardInterrupt
Child process 68807 (/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/recover.lock.QC_hg38__K562_RNA-seq_60_bamQC.tsv

### Pipeline failed at:  (06-11 18:55:15) elapsed: 0.0 _TIME_

Total time: 1:39:58
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_60 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_60pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba25-33c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/
*  Pipeline started at:   (06-11 19:08:42) elapsed: 2.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_60pctRNA.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_60`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_60pctRNA.fastq.gz

> `File_mb`	5477.44	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 19:08:43) elapsed: 2.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/K562_RNA-seq_60.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/K562_RNA-seq_60.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/K562_RNA-seq_60_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/K562_RNA-seq_60.fastq.gz']

### FASTQ processing:  (06-11 19:08:43) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-11 19:08:44) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt/K562_RNA-seq_60_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_60_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_60_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-11 19:08:44) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 19:08:44) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/prealignments/K562_RNA-seq_60_human_rDNA_unmap.fq

### Map to genome (06-11 19:08:44) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam`  

### Compress all unmapped read files (06-11 19:08:44) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/prealignments/K562_RNA-seq_60_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-11 19:08:44) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam.bai`  
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/lock.QC_hg38__K562_RNA-seq_60_bamQC.tsv
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_bamQC.tsv` (5131)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/tmp_K562_RNA-seq_60_sort_mxgcxuey'
Processing with 12 cores...
Discarding 98 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 97 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:01:10. Running peak memory: 1.352GB.  
  PID: 5131;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.352GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_bamQC.tsv`

> `NRF`	0.69	PEPPRO	_RES_

> `PBC1`	0.89	PEPPRO	_RES_

> `PBC2`	14.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_unmap.bam`  

> `samtools view -b -@ 12 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_unmap.bam` (5234)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 1.352GB.  
  PID: 5234;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_temp.bam`

> `Unmapped_reads`	7219211	PEPPRO	_RES_

### Split BAM by strand (06-11 19:10:36) elapsed: 112.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam` (5595)
<pre>
</pre>
Command completed. Elapsed time: 0:03:43. Running peak memory: 1.352GB.  
  PID: 5595;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_minus.bam` (5800)
<pre>
</pre>
Command completed. Elapsed time: 0:03:36. Running peak memory: 1.352GB.  
  PID: 5800;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-11 19:17:55) elapsed: 439.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (6201)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.352GB.  
  PID: 6201;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_plus_TssEnrichment.txt` (6202)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 1.352GB.  
  PID: 6202;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.524GB


> `TSS_coding_score`	17.4	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_minus_TssEnrichment.txt` (6238)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 1.352GB.  
  PID: 6238;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.535GB


> `TSS_non-coding_score`	4.2	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_minus_TssEnrichment.txt` (6271)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 1.352GB.  
  PID: 6271;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_60_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_60_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt` (6295,6296,6297,6298)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.352GB.  
  PID: 6295;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 6297;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 6296;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 6298;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_keep.txt` (6301)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.352GB.  
  PID: 6301;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-11 19:18:18) elapsed: 23.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_ensembl_tss.bed` (6303,6304)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 1.352GB.  
  PID: 6303;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 6304;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_ensembl_gene_body.bed` (6307,6308)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.352GB.  
  PID: 6307;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 6308;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_TSS_density.bed` (6310,6311,6312,6313)
<pre>
</pre>
Command completed. Elapsed time: 0:01:11. Running peak memory: 1.352GB.  
  PID: 6311;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 6313;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 6310;	Command: bedtools;	Return code: 0;	Memory used: 0.046GB  
  PID: 6312;	Command: sort;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_gene_body_density.bed` (6378,6379,6380)
<pre>
</pre>
Command completed. Elapsed time: 0:01:37. Running peak memory: 1.352GB.  
  PID: 6380;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 6378;	Command: bedtools;	Return code: 0;	Memory used: 0.164GB  
  PID: 6379;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/tmpg_a6oc59` (6749,6750,6751)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.352GB.  
  PID: 6749;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 6751;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 6750;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/tmpg_a6oc59 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-11 19:21:08) elapsed: 170.0 _TIME_

Total time: 0:12:28
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/tmpg_a6oc59 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_60 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_60pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aw29-25b
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/
*  Pipeline started at:   (06-14 21:11:48) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_60pctRNA.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_60`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_60pctRNA.fastq.gz

> `File_mb`	5477.44	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-14 21:11:48) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/K562_RNA-seq_60.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/K562_RNA-seq_60.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/K562_RNA-seq_60_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/K562_RNA-seq_60.fastq.gz']

### FASTQ processing:  (06-14 21:11:48) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-14 21:11:48) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/cutadapt/K562_RNA-seq_60_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_60_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_60_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-14 21:11:48) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:11:48) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/prealignments/K562_RNA-seq_60_human_rDNA_unmap.fq

### Map to genome (06-14 21:11:48) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam`  

### Compress all unmapped read files (06-14 21:11:48) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/prealignments/K562_RNA-seq_60_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:11:48) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_unmap.bam`  

### Split BAM by strand (06-14 21:11:48) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_minus.bam`  

### Calculate TSS enrichment (06-14 21:11:48) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/K562_RNA-seq_60_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_60_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:11:48) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/tmpo2ptjhlp` (379316,379317,379318)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.006GB.  
  PID: 379316;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 379318;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 379317;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/tmpo2ptjhlp | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0390275) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/tmpo2ptjhlp > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_pause_index.bed` (379324)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.006GB.  
  PID: 379324;	Command: awk;	Return code: 0;	Memory used: 0.004GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	8.19	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_pause_index.bed` (379329)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.319GB.  
  PID: 379329;	Command: Rscript;	Return code: 0;	Memory used: 0.319GB

> `Pause index`	QC_hg38/K562_RNA-seq_60_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_60_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_pause_index.bed` (379350)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.319GB.  
  PID: 379350;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:11:54) elapsed: 5.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam`
48343348 18391165

> `Plus_FRiP`	0.38	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_minus.bam`
48343348 18206742

> `Minus_FRiP`	0.38	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/signal_hg38/K562_RNA-seq_60_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_gene_sort.bed` (387365,387369)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.319GB.  
  PID: 387365;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 387369;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/signal_hg38/K562_RNA-seq_60_gene_coverage.bed` (387779)
<pre>
</pre>
Command completed. Elapsed time: 0:01:36. Running peak memory: 0.319GB.  
  PID: 387779;	Command: bedtools;	Return code: 0;	Memory used: 0.175GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/hg38_annotations.bed.gz` (407094)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.319GB.  
  PID: 407094;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/hg38_annotations.bed` (407095)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.319GB.  
  PID: 407095;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:14:51) elapsed: 177.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/raw/hg38_annotations.bed` (407103)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.319GB.  
  PID: 407103;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Enhancer_sort.bed` (407105,407106,407107,407108)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.319GB.  
  PID: 407105;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 407106;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 407108;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB  
  PID: 407107;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Enhancer_plus_coverage.bed` (407111)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 0.319GB.  
  PID: 407111;	Command: bedtools;	Return code: 0;	Memory used: 0.025GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Enhancer_minus_coverage.bed` (407349)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 0.319GB.  
  PID: 407349;	Command: bedtools;	Return code: 0;	Memory used: 0.032GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter_sort.bed` (407397,407398,407399,407400)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.319GB.  
  PID: 407397;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 407398;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 407400;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 407399;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Promoter_plus_coverage.bed` (407402)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 0.319GB.  
  PID: 407402;	Command: bedtools;	Return code: 0;	Memory used: 0.124GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Promoter_minus_coverage.bed` (407456)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 0.319GB.  
  PID: 407456;	Command: bedtools;	Return code: 0;	Memory used: 0.065GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter_Flanking_Region"` (410886)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.319GB.  
  PID: 410886;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter_Flanking_Region_sort.bed` (410887,410888,410889,410890)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.319GB.  
  PID: 410887;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 410889;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 410888;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 410890;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Promoter_Flanking_Region_plus_coverage.bed` (411111)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 0.319GB.  
  PID: 411111;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Promoter_Flanking_Region_minus_coverage.bed` (413063)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 0.319GB.  
  PID: 413063;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/5_UTR"` (421839)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.319GB.  
  PID: 421839;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/5_UTR_sort.bed` (421896,421897,421898,421900)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.319GB.  
  PID: 421896;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 421897;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 421900;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 421898;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_5_UTR_plus_coverage.bed` (422372)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 0.319GB.  
  PID: 422372;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_5_UTR_minus_coverage.bed` (429656)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 0.319GB.  
  PID: 429656;	Command: bedtools;	Return code: 0;	Memory used: 0.026GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/3_UTR"` (429712)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.319GB.  
  PID: 429712;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/3_UTR_sort.bed` (429713,429714,429715,429716)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.319GB.  
  PID: 429713;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 429714;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 429716;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB  
  PID: 429715;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_3_UTR_plus_coverage.bed` (429719)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 0.319GB.  
  PID: 429719;	Command: bedtools;	Return code: 0;	Memory used: 0.04GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_3_UTR_minus_coverage.bed` (430127)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 0.319GB.  
  PID: 430127;	Command: bedtools;	Return code: 0;	Memory used: 0.05GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Exon_sort.bed` (430217,430218,430219,430220)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.319GB.  
  PID: 430217;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 430218;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 430220;	Command: bedtools;	Return code: 0;	Memory used: 0.16GB  
  PID: 430219;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Exon_plus_coverage.bed` (430224)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 0.319GB.  
  PID: 430224;	Command: bedtools;	Return code: 0;	Memory used: 0.131GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Exon_minus_coverage.bed` (430308)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 0.319GB.  
  PID: 430308;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Intron_sort.bed` (438319,438331,438339,438343)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.319GB.  
  PID: 438319;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 438339;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 438331;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 438343;	Command: bedtools;	Return code: 0;	Memory used: 0.078GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Intron_plus_coverage.bed` (439176)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 0.319GB.  
  PID: 439176;	Command: bedtools;	Return code: 0;	Memory used: 0.116GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Intron_minus_coverage.bed` (447587)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 0.319GB.  
  PID: 447587;	Command: bedtools;	Return code: 0;	Memory used: 0.06GB


### Plot cFRiF/FRiF (06-14 21:23:31) elapsed: 520.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_RNA-seq_60 -z 3099922541 -n 23047249 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Intron_plus_coverage.bed` (453190)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 0.461GB.  
  PID: 453190;	Command: Rscript;	Return code: 0;	Memory used: 0.461GB

> `cFRiF`	QC_hg38/K562_RNA-seq_60_cFRiF.pdf	cFRiF	QC_hg38/K562_RNA-seq_60_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_RNA-seq_60 -z 3099922541 -n 23047249 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_Intron_plus_coverage.bed` (453286)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 0.461GB.  
  PID: 453286;	Command: Rscript;	Return code: 0;	Memory used: 0.461GB

> `FRiF`	QC_hg38/K562_RNA-seq_60_FRiF.pdf	FRiF	QC_hg38/K562_RNA-seq_60_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 21:24:38) elapsed: 67.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_exons_sort.bed` (453340,453342)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.461GB.  
  PID: 453342;	Command: bedtools;	Return code: 0;	Memory used: 0.094GB  
  PID: 453340;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_introns_sort.bed` (453348,453349,453350)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.461GB.  
  PID: 453348;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 453350;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 453349;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_exons_coverage.bed` (453437)
<pre>
</pre>
Command completed. Elapsed time: 0:01:20. Running peak memory: 0.461GB.  
  PID: 453437;	Command: bedtools;	Return code: 0;	Memory used: 0.144GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_introns_coverage.bed` (17387)
<pre>
</pre>
Command completed. Elapsed time: 0:01:20. Running peak memory: 0.461GB.  
  PID: 17387;	Command: bedtools;	Return code: 0;	Memory used: 0.132GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/48.343348)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_exons_rpkm.bed` (18107,18108,18109)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.461GB.  
  PID: 18107;	Command: awk;	Return code: 0;	Memory used: 0.009GB  
  PID: 18109;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 18108;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/48.343348)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_introns_rpkm.bed` (18112,18113,18114)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.461GB.  
  PID: 18112;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 18114;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 18113;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_exon_intron_ratios.bed` (18117,18118,18119)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.461GB.  
  PID: 18117;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 18119;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 18118;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	4.92	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_exon_intron_ratios.bed --annotate` (18126)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.461GB.  
  PID: 18126;	Command: Rscript;	Return code: 0;	Memory used: 0.317GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_60_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_60_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/QC_hg38/K562_RNA-seq_60_exon_intron_ratios.bed` (18159)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.461GB.  
  PID: 18159;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (06-14 21:27:36) elapsed: 178.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/signal_hg38/K562_RNA-seq_60_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/signal_hg38/K562_RNA-seq_60_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam` (18167)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 0.461GB.  
  PID: 18167;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/signal_hg38/K562_RNA-seq_60_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/signal_hg38/K562_RNA-seq_60_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 48343348.0` (18245)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_60_plus_cuttrace_zyw3xwxr'
Processing with 4 cores...
stdin is empty of data
stdin is empty of data
Discarding 108 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 87 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 87 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/signal_hg38/K562_RNA-seq_60_plus_exact_body_0-mer.bw'
Merging 87 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/signal_hg38/K562_RNA-seq_60_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:06. Running peak memory: 2.591GB.  
  PID: 18245;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.591GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/signal_hg38/K562_RNA-seq_60_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/signal_hg38/K562_RNA-seq_60_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_minus.bam` (43448)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 2.591GB.  
  PID: 43448;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/signal_hg38/K562_RNA-seq_60_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/signal_hg38/K562_RNA-seq_60_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 48343348.0` (43471)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/aligned_hg38/K562_RNA-seq_60_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_60_minus_cuttrace_a6cvchqb'
Processing with 4 cores...
stdin is empty of data
Discarding 112 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1']
Keeping 83 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/signal_hg38/K562_RNA-seq_60_minus_exact_body_0-mer.bw'
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_60/signal_hg38/K562_RNA-seq_60_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:50. Running peak memory: 2.591GB.  
  PID: 43471;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.504GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:32:30
*  Total elapsed time (all runs):  2:38:26
*         Peak memory (this run):  2.5912 GB
*        Pipeline completed time: 2020-06-14 21:44:18
