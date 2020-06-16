### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_60 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_60pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 24 -M 24000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aw29-23a
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/
*  Pipeline started at:   (06-15 07:11:05) elapsed: 2.0 _TIME_

### Version log:

*       Python version:  3.6.6
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.9.8
*        Pipeline hash:  887a9a85408962dc3ca070dc954e59a5d4d73a10
*      Pipeline branch:  * master
*        Pipeline date:  2020-06-09 11:36:49 -0400
*        Pipeline diff:  40 files changed, 16413 insertions(+), 3702 deletions(-)

### Arguments passed to pipeline:

*           `TSS_name`:  `None`
*            `adapter`:  `cutadapt`
*          `anno_name`:  `None`
*         `complexity`:  `False`
*        `config_file`:  `peppro.yaml`
*              `cores`:  `24`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_60pct.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `24000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `K562_PRO-seq_60`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_60pct.fastq.gz

> `File_mb`	21820.91	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:11:05) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/raw/K562_PRO-seq_60.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_PRO_60pct.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/raw/K562_PRO-seq_60.fastq.gz` (12986)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 12986;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/raw/K562_PRO-seq_60.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastq/K562_PRO-seq_60_R1.fastq`  

> `pigz -f -p 24 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/raw/K562_PRO-seq_60.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastq/K562_PRO-seq_60_R1.fastq` (12987)
<pre>
</pre>
Command completed. Elapsed time: 0:10:12. Running peak memory: 0.003GB.  
  PID: 12987;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	297971266	PEPPRO	_RES_

> `Fastq_reads`	297971266	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/raw/K562_PRO-seq_60.fastq.gz']

### FASTQ processing:  (06-15 07:36:03) elapsed: 1498.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastq/K562_PRO-seq_60_R1_processed.fastq`  

> `(cutadapt -j 24 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastq/K562_PRO-seq_60_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastq/K562_PRO-seq_60_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastq/K562_PRO-seq_60_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/cutadapt/K562_PRO-seq_60_R1_cutadapt.txt` (15755)
<pre>
</pre>
Command completed. Elapsed time: 0:10:57. Running peak memory: 9.924GB.  
  PID: 15755;	Command: cutadapt;	Return code: 0;	Memory used: 9.924GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastq/K562_PRO-seq_60_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastq/K562_PRO-seq_60_R1_processed.fastq` (16816,16817)
<pre>
</pre>
Command completed. Elapsed time: 0:07:12. Running peak memory: 9.924GB.  
  PID: 16816;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 16817;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/cutadapt/K562_PRO-seq_60_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	253852065.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/cutadapt/K562_PRO-seq_60_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/cutadapt/K562_PRO-seq_60_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastq/K562_PRO-seq_60_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	7461091.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	2.504	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastqc/K562_PRO-seq_60_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (17874)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 17874;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	290510175	PEPPRO	_RES_

> `Trim_loss_rate`	2.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastq/K562_PRO-seq_60_R1_processed.fastq` (17939)
<pre>
Started analysis of K562_PRO-seq_60_R1_processed.fastq
Approx 5% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 10% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 15% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 20% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 25% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 30% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 35% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 40% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 45% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 50% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 55% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 60% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 65% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 70% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 75% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 80% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 85% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 90% complete for K562_PRO-seq_60_R1_processed.fastq
Approx 95% complete for K562_PRO-seq_60_R1_processed.fastq
Analysis complete for K562_PRO-seq_60_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:15:57. Running peak memory: 9.924GB.  
  PID: 17939;	Command: fastqc;	Return code: 0;	Memory used: 0.261GB

> `FastQC report r1`	fastqc/K562_PRO-seq_60_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastq/processed_R1.flag` (19680)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 19680;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-15 08:15:42) elapsed: 2379.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/cutadapt/K562_PRO-seq_60_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/cutadapt/K562_PRO-seq_60_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/cutadapt` (19686)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 9.924GB.  
  PID: 19686;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_60_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_60_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/cutadapt/K562_PRO-seq_60_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 08:15:47) elapsed: 5.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/cutadapt/K562_PRO-seq_60_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/cutadapt/K562_PRO-seq_60_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/cutadapt/K562_PRO-seq_60_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/cutadapt/K562_PRO-seq_60_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/cutadapt/K562_PRO-seq_60_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2313	PEPPRO	_RES_

### Prealignments (06-15 08:15:47) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 08:15:47) elapsed: 0.0 _TIME_


> `(bowtie2 -p 24 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_PRO-seq_60 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/fastq/K562_PRO-seq_60_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/prealignments/K562_PRO-seq_60_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
290510175 reads; of these:
  290510175 (100.00%) were unpaired; of these:
    263700144 (90.77%) aligned 0 times
    26810031 (9.23%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.23% overall alignment rate

> `Aligned_reads_human_rDNA`	26810031.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.23	PEPPRO	_RES_

### Map to genome (06-15 08:51:53) elapsed: 2166.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam`  

> `bowtie2 -p 24 --very-sensitive --rg-id K562_PRO-seq_60 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/prealignments/K562_PRO-seq_60_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/tmpw8xzl21o -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_temp.bam` (23413,23414,23415)
<pre>
263700144 reads; of these:
  263700144 (100.00%) were unpaired; of these:
    3235155 (1.23%) aligned 0 times
    196904698 (74.67%) aligned exactly 1 time
    63560291 (24.10%) aligned >1 times
98.77% overall alignment rate
[bam_sort_core] merging from 84 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 2:06:28. Running peak memory: 9.924GB.  
  PID: 23413;	Command: bowtie2;	Return code: 0;	Memory used: 4.009GB  
  PID: 23414;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 23415;	Command: samtools;	Return code: 0;	Memory used: 0.901GB


> `samtools view -q 10 -b -@ 24 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam` (80374)
<pre>
</pre>
Command completed. Elapsed time: 0:06:18. Running peak memory: 9.924GB.  
  PID: 80374;	Command: samtools;	Return code: 0;	Memory used: 0.029GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	260464989	PEPPRO	_RES_

> `QC_filtered_reads`	28049270	PEPPRO	_RES_

> `Aligned_reads`	232415719	PEPPRO	_RES_

> `Alignment_rate`	80.0	PEPPRO	_RES_

> `Total_efficiency`	78.0	PEPPRO	_RES_

> `Read_depth`	18.77	PEPPRO	_RES_

### Compress all unmapped read files (06-15 11:46:46) elapsed: 10492.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/prealignments/K562_PRO-seq_60_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 24 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/prealignments/K562_PRO-seq_60_human_rDNA_unmap.fq` (95707)
<pre>
</pre>
Command completed. Elapsed time: 0:04:56. Running peak memory: 9.924GB.  
  PID: 95707;	Command: pigz;	Return code: 0;	Memory used: 0.019GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_temp.bam` (96244)
<pre>
</pre>
Command completed. Elapsed time: 0:04:01. Running peak memory: 9.924GB.  
  PID: 96244;	Command: samtools;	Return code: 0;	Memory used: 0.025GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	5490660	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam` (96651)
<pre>
</pre>
Command completed. Elapsed time: 0:03:18. Running peak memory: 9.924GB.  
  PID: 96651;	Command: samtools;	Return code: 0;	Memory used: 0.021GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/chr_sizes.bed` (96824,96825,96826,96827)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 96826;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 96824;	Command: samtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 96827;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 96825;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/chr_sizes.bed -b -@ 24 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_noMT.bam` (96829)
<pre>
</pre>
Command completed. Elapsed time: 0:03:33. Running peak memory: 9.924GB.  
  PID: 96829;	Command: samtools;	Return code: 0;	Memory used: 0.03GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam` (97305)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 9.924GB.  
  PID: 97305;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam` (97308)
<pre>
</pre>
Command completed. Elapsed time: 0:03:21. Running peak memory: 9.924GB.  
  PID: 97308;	Command: samtools;	Return code: 0;	Memory used: 0.02GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 12:16:14) elapsed: 1768.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam -c 24 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_bamQC.tsv` (98683)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/tmp_K562_PRO-seq_60_sort_ylp8f_g6'
Processing with 24 cores...
Discarding 86 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 109 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:03:53. Running peak memory: 9.924GB.  
  PID: 98683;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 7.598GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_bamQC.tsv`

> `NRF`	0.42	PEPPRO	_RES_

> `PBC1`	0.66	PEPPRO	_RES_

> `PBC2`	4.61	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_unmap.bam`  

> `samtools view -b -@ 24 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_unmap.bam` (99172)
<pre>
</pre>
Command completed. Elapsed time: 0:00:47. Running peak memory: 9.924GB.  
  PID: 99172;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools view -c -f 4 -@ 24 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_temp.bam`

> `Unmapped_reads`	3235155	PEPPRO	_RES_

### Split BAM by strand (06-15 12:21:34) elapsed: 321.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_plus.bam` (99310)
<pre>
</pre>
Command completed. Elapsed time: 0:16:25. Running peak memory: 9.924GB.  
  PID: 99310;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_minus.bam` (100828)
<pre>
</pre>
Command completed. Elapsed time: 0:15:40. Running peak memory: 9.924GB.  
  PID: 100828;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 12:53:39) elapsed: 1925.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (102319)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 102319;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/plus_TSS.tsv -p ends -c 24 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_plus_TssEnrichment.txt` (102326)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 9.924GB.  
  PID: 102326;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.555GB


> `TSS_coding_score`	13.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/minus_TSS.tsv -p ends -c 24 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_minus_TssEnrichment.txt` (102397)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 9.924GB.  
  PID: 102397;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.718GB


> `TSS_non-coding_score`	5.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_minus_TssEnrichment.txt` (102468)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 9.924GB.  
  PID: 102468;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `TSS enrichment`	QC_hg38/K562_PRO-seq_60_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_60_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt` (102496,102497,102498,102499)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 102496;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 102498;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 102497;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 102499;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_keep.txt` (102501)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 102501;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 12:54:29) elapsed: 50.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/hg38_ensembl_tss.bed` (102503,102504)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 9.924GB.  
  PID: 102503;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 102504;	Command: bedtools;	Return code: 0;	Memory used: 0.096GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/hg38_ensembl_gene_body.bed` (102508,102509)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 102508;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 102509;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_TSS_density.bed` (102511,102512,102513,102514)
<pre>
</pre>
Command completed. Elapsed time: 0:04:56. Running peak memory: 9.924GB.  
  PID: 102511;	Command: bedtools;	Return code: 0;	Memory used: 0.064GB  
  PID: 102513;	Command: sort;	Return code: 0;	Memory used: 0.012GB  
  PID: 102512;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 102514;	Command: sort;	Return code: 0;	Memory used: 0.015GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_gene_body_density.bed` (102964,102965,102966)
<pre>
</pre>
Command completed. Elapsed time: 0:08:08. Running peak memory: 9.924GB.  
  PID: 102965;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 102964;	Command: bedtools;	Return code: 0;	Memory used: 0.623GB  
  PID: 102966;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/tmp6hgx34r6` (103875,103882,103883)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 9.924GB.  
  PID: 103875;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 103883;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 103882;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/tmp6hgx34r6 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.137871) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/tmp6hgx34r6 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_pause_index.bed` (103889)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 103889;	Command: awk;	Return code: 0;	Memory used: 0.004GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.88	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_pause_index.bed` (103894)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 9.924GB.  
  PID: 103894;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Pause index`	QC_hg38/K562_PRO-seq_60_pause_index.pdf	Pause index	QC_hg38/K562_PRO-seq_60_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_pause_index.bed.gz`  

> `pigz -f -p 24 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_pause_index.bed` (103921)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 103921;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 13:07:41) elapsed: 792.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_plus.bam`
232415719 82895263

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_minus.bam`
232415719 78848428

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/signal_hg38/K562_PRO-seq_60_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/hg38_gene_sort.bed` (104510,104511)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 9.924GB.  
  PID: 104510;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 104511;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/signal_hg38/K562_PRO-seq_60_gene_coverage.bed` (104514)
<pre>
</pre>
Command completed. Elapsed time: 0:08:07. Running peak memory: 9.924GB.  
  PID: 104514;	Command: bedtools;	Return code: 0;	Memory used: 0.621GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/raw/hg38_annotations.bed.gz` (105370)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 105370;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 24 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/raw/hg38_annotations.bed` (105375)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 105375;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 13:22:19) elapsed: 878.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/raw/hg38_annotations.bed` (105386)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 9.924GB.  
  PID: 105386;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Enhancer_sort.bed` (105388,105389,105390,105391)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 9.924GB.  
  PID: 105388;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 105389;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 105391;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB  
  PID: 105390;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Enhancer_plus_coverage.bed` (105394)
<pre>
</pre>
Command completed. Elapsed time: 0:02:26. Running peak memory: 9.924GB.  
  PID: 105394;	Command: bedtools;	Return code: 0;	Memory used: 0.025GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Enhancer_minus_coverage.bed` (105524)
<pre>
</pre>
Command completed. Elapsed time: 0:02:23. Running peak memory: 9.924GB.  
  PID: 105524;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter_sort.bed` (105839,105840,105841,105842)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 105839;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 105840;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 105842;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 105841;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Promoter_plus_coverage.bed` (105844)
<pre>
</pre>
Command completed. Elapsed time: 0:02:42. Running peak memory: 9.924GB.  
  PID: 105844;	Command: bedtools;	Return code: 0;	Memory used: 0.336GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Promoter_minus_coverage.bed` (106162)
<pre>
</pre>
Command completed. Elapsed time: 0:02:42. Running peak memory: 9.924GB.  
  PID: 106162;	Command: bedtools;	Return code: 0;	Memory used: 0.286GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter_Flanking_Region"` (106547)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 106547;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter_Flanking_Region_sort.bed` (106549,106550,106551,106552)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 9.924GB.  
  PID: 106549;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 106551;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 106550;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 106552;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Promoter_Flanking_Region_plus_coverage.bed` (106554)
<pre>
</pre>
Command completed. Elapsed time: 0:02:33. Running peak memory: 9.924GB.  
  PID: 106554;	Command: bedtools;	Return code: 0;	Memory used: 0.062GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Promoter_Flanking_Region_minus_coverage.bed` (106946)
<pre>
</pre>
Command completed. Elapsed time: 0:02:30. Running peak memory: 9.924GB.  
  PID: 106946;	Command: bedtools;	Return code: 0;	Memory used: 0.102GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/5_UTR"` (107104)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 107104;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/5_UTR_sort.bed` (107105,107106,107107,107108)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 9.924GB.  
  PID: 107105;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 107106;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 107108;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 107107;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_5_UTR_plus_coverage.bed` (107111)
<pre>
</pre>
Command completed. Elapsed time: 0:02:29. Running peak memory: 9.924GB.  
  PID: 107111;	Command: bedtools;	Return code: 0;	Memory used: 0.048GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_5_UTR_minus_coverage.bed` (107558)
<pre>
</pre>
Command completed. Elapsed time: 0:02:23. Running peak memory: 9.924GB.  
  PID: 107558;	Command: bedtools;	Return code: 0;	Memory used: 0.032GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/3_UTR"` (107689)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 107689;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/3_UTR_sort.bed` (107690,107691,107692,107693)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 9.924GB.  
  PID: 107690;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 107691;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 107693;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB  
  PID: 107692;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_3_UTR_plus_coverage.bed` (107695)
<pre>
</pre>
Command completed. Elapsed time: 0:02:30. Running peak memory: 9.924GB.  
  PID: 107695;	Command: bedtools;	Return code: 0;	Memory used: 0.065GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_3_UTR_minus_coverage.bed` (108018)
<pre>
</pre>
Command completed. Elapsed time: 0:02:26. Running peak memory: 9.924GB.  
  PID: 108018;	Command: bedtools;	Return code: 0;	Memory used: 0.1GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Exon_sort.bed` (108153,108154,108155,108156)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 9.924GB.  
  PID: 108153;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 108154;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 108156;	Command: bedtools;	Return code: 0;	Memory used: 0.162GB  
  PID: 108155;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Exon_plus_coverage.bed` (108161)
<pre>
</pre>
Command completed. Elapsed time: 0:02:41. Running peak memory: 9.924GB.  
  PID: 108161;	Command: bedtools;	Return code: 0;	Memory used: 0.283GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Exon_minus_coverage.bed` (108606)
<pre>
</pre>
Command completed. Elapsed time: 0:02:35. Running peak memory: 9.924GB.  
  PID: 108606;	Command: bedtools;	Return code: 0;	Memory used: 0.114GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Intron_sort.bed` (108743,108744,108745,108746)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 9.924GB.  
  PID: 108743;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 108745;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 108744;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 108746;	Command: bedtools;	Return code: 0;	Memory used: 0.083GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Intron_plus_coverage.bed` (108749)
<pre>
</pre>
Command completed. Elapsed time: 0:03:09. Running peak memory: 9.924GB.  
  PID: 108749;	Command: bedtools;	Return code: 0;	Memory used: 0.277GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Intron_minus_coverage.bed` (109110)
<pre>
</pre>
Command completed. Elapsed time: 0:03:11. Running peak memory: 9.924GB.  
  PID: 109110;	Command: bedtools;	Return code: 0;	Memory used: 0.266GB


### Plot cFRiF/FRiF (06-15 13:59:09) elapsed: 2210.0 _TIME_


> `samtools view -@ 24 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_60 -z 3099922541 -n 115653622 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Intron_plus_coverage.bed` (109323)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 9.924GB.  
  PID: 109323;	Command: Rscript;	Return code: 0;	Memory used: 0.451GB

> `cFRiF`	QC_hg38/K562_PRO-seq_60_cFRiF.pdf	cFRiF	QC_hg38/K562_PRO-seq_60_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_60 -z 3099922541 -n 115653622 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_Intron_plus_coverage.bed` (109685)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 9.924GB.  
  PID: 109685;	Command: Rscript;	Return code: 0;	Memory used: 0.457GB

> `FRiF`	QC_hg38/K562_PRO-seq_60_FRiF.pdf	FRiF	QC_hg38/K562_PRO-seq_60_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 14:00:35) elapsed: 86.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/hg38_exons_sort.bed` (109722,109723)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 9.924GB.  
  PID: 109723;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB  
  PID: 109722;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/hg38_introns_sort.bed` (109730,109731,109732)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 9.924GB.  
  PID: 109730;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 109732;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 109731;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_exons_coverage.bed` (109738)
<pre>
</pre>
Command completed. Elapsed time: 0:05:18. Running peak memory: 9.924GB.  
  PID: 109738;	Command: bedtools;	Return code: 0;	Memory used: 0.179GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_introns_coverage.bed` (110256)
<pre>
</pre>
Command completed. Elapsed time: 0:06:44. Running peak memory: 9.924GB.  
  PID: 110256;	Command: bedtools;	Return code: 0;	Memory used: 0.289GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/232.415719)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_exons_rpkm.bed` (110924,110935,110936)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 9.924GB.  
  PID: 110924;	Command: awk;	Return code: 0;	Memory used: 0.011GB  
  PID: 110936;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 110935;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/232.415719)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_introns_rpkm.bed` (110938,110939,110940)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 9.924GB.  
  PID: 110938;	Command: awk;	Return code: 0;	Memory used: 0.025GB  
  PID: 110940;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 110939;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_exon_intron_ratios.bed` (110943,110944,110945)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 110943;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 110945;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 110944;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_exon_intron_ratios.bed --annotate` (110951)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 9.924GB.  
  PID: 110951;	Command: Rscript;	Return code: 0;	Memory used: 0.11GB

> `mRNA contamination`	QC_hg38/K562_PRO-seq_60_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_PRO-seq_60_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_exon_intron_ratios.bed.gz`  

> `pigz -f -p 24 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/QC_hg38/K562_PRO-seq_60_exon_intron_ratios.bed` (110982)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 9.924GB.  
  PID: 110982;	Command: pigz;	Return code: 0;	Memory used: 0.001GB


### Produce bigWig files (06-15 14:13:05) elapsed: 750.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/signal_hg38/K562_PRO-seq_60_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/signal_hg38/K562_PRO-seq_60_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_plus.bam` (110992)
<pre>
</pre>
Command completed. Elapsed time: 0:01:42. Running peak memory: 9.924GB.  
  PID: 110992;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/signal_hg38/K562_PRO-seq_60_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/signal_hg38/K562_PRO-seq_60_plus_smooth_body_0-mer.bw -p 16 --variable-step --tail-edge --scale 232415719.0` (111082)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_plus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_60_plus_cuttrace_hf8208zn'
Processing with 8 cores...
stdin is empty of data
psutil.NoSuchProcess process no longer exists (pid=111684)
Warning: couldn't add memory use for process: 111082
Discarding 92 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1']
Keeping 103 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 103 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/signal_hg38/K562_PRO-seq_60_plus_exact_body_0-mer.bw'
Merging 103 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/signal_hg38/K562_PRO-seq_60_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:22. Running peak memory: 9.924GB.  
  PID: 111082;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.505GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/signal_hg38/K562_PRO-seq_60_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/signal_hg38/K562_PRO-seq_60_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_minus.bam` (112997)
<pre>
</pre>
Command completed. Elapsed time: 0:01:39. Running peak memory: 9.924GB.  
  PID: 112997;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/signal_hg38/K562_PRO-seq_60_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/signal_hg38/K562_PRO-seq_60_minus_smooth_body_0-mer.bw -p 16 --variable-step --tail-edge --scale 232415719.0` (113234)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/aligned_hg38/K562_PRO-seq_60_minus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_60_minus_cuttrace__hslzqh6'
Processing with 8 cores...
stdin is empty of data
stdin is empty of data
Discarding 97 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270723v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 98 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 98 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/signal_hg38/K562_PRO-seq_60_minus_exact_body_0-mer.bw'
Merging 98 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_60/signal_hg38/K562_PRO-seq_60_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:01. Running peak memory: 9.924GB.  
  PID: 113234;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.729GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  7:21:47
*  Total elapsed time (all runs):  9:06:47
*         Peak memory (this run):  9.9242 GB
*        Pipeline completed time: 2020-06-15 14:32:49
