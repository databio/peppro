### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_20 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_20pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aw29-25a
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/
*  Pipeline started at:   (06-15 06:58:48) elapsed: 57.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_20pct.fastq.gz']`
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
*        `sample_name`:  `K562_PRO-seq_20`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_20pct.fastq.gz

> `File_mb`	7335.34	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 06:58:49) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/raw/K562_PRO-seq_20.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_PRO_20pct.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/raw/K562_PRO-seq_20.fastq.gz` (261342)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.002GB.  
  PID: 261342;	Command: ln;	Return code: 0;	Memory used: 0.002GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/raw/K562_PRO-seq_20.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastq/K562_PRO-seq_20_R1.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/raw/K562_PRO-seq_20.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastq/K562_PRO-seq_20_R1.fastq` (261343)
<pre>
</pre>
Command completed. Elapsed time: 0:02:41. Running peak memory: 0.003GB.  
  PID: 261343;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	99310774	PEPPRO	_RES_

> `Fastq_reads`	99310774	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/raw/K562_PRO-seq_20.fastq.gz']

### FASTQ processing:  (06-15 07:03:33) elapsed: 284.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastq/K562_PRO-seq_20_R1_processed.fastq`  

> `(cutadapt -j 12 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastq/K562_PRO-seq_20_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastq/K562_PRO-seq_20_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastq/K562_PRO-seq_20_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/cutadapt/K562_PRO-seq_20_R1_cutadapt.txt` (261906)
<pre>
</pre>
Command completed. Elapsed time: 0:02:36. Running peak memory: 4.632GB.  
  PID: 261906;	Command: cutadapt;	Return code: 0;	Memory used: 4.632GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastq/K562_PRO-seq_20_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastq/K562_PRO-seq_20_R1_processed.fastq` (262260,262261)
<pre>
</pre>
Command completed. Elapsed time: 0:01:58. Running peak memory: 4.632GB.  
  PID: 262260;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 262261;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/cutadapt/K562_PRO-seq_20_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	84602654.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/cutadapt/K562_PRO-seq_20_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/cutadapt/K562_PRO-seq_20_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastq/K562_PRO-seq_20_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	2486025.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	2.5033	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastqc/K562_PRO-seq_20_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (262414)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 262414;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	96824749	PEPPRO	_RES_

> `Trim_loss_rate`	2.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastq/K562_PRO-seq_20_R1_processed.fastq` (262435)
<pre>
Started analysis of K562_PRO-seq_20_R1_processed.fastq
Approx 5% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 10% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 15% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 20% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 25% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 30% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 35% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 40% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 45% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 50% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 55% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 60% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 65% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 70% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 75% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 80% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 85% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 90% complete for K562_PRO-seq_20_R1_processed.fastq
Approx 95% complete for K562_PRO-seq_20_R1_processed.fastq
Analysis complete for K562_PRO-seq_20_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:05:41. Running peak memory: 4.632GB.  
  PID: 262435;	Command: fastqc;	Return code: 0;	Memory used: 0.216GB

> `FastQC report r1`	fastqc/K562_PRO-seq_20_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastq/processed_R1.flag` (262986)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 262986;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-15 07:14:50) elapsed: 677.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/cutadapt/K562_PRO-seq_20_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/cutadapt/K562_PRO-seq_20_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/cutadapt` (262992)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 4.632GB.  
  PID: 262992;	Command: Rscript;	Return code: 0;	Memory used: 0.123GB

> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_20_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_20_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/cutadapt/K562_PRO-seq_20_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:14:58) elapsed: 8.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/cutadapt/K562_PRO-seq_20_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/cutadapt/K562_PRO-seq_20_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/cutadapt/K562_PRO-seq_20_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/cutadapt/K562_PRO-seq_20_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/cutadapt/K562_PRO-seq_20_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2311	PEPPRO	_RES_

### Prealignments (06-15 07:14:58) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:14:58) elapsed: 0.0 _TIME_


> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_PRO-seq_20 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/fastq/K562_PRO-seq_20_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/prealignments/K562_PRO-seq_20_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
96824749 reads; of these:
  96824749 (100.00%) were unpaired; of these:
    87890105 (90.77%) aligned 0 times
    8934644 (9.23%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.23% overall alignment rate

> `Aligned_reads_human_rDNA`	8934644.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.23	PEPPRO	_RES_

### Map to genome (06-15 07:30:24) elapsed: 925.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam`  

> `bowtie2 -p 12 --very-sensitive --rg-id K562_PRO-seq_20 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/prealignments/K562_PRO-seq_20_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/tmprvjzknrw -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_temp.bam` (264826,264831,264832)
<pre>
87890105 reads; of these:
  87890105 (100.00%) were unpaired; of these:
    1079115 (1.23%) aligned 0 times
    65634849 (74.68%) aligned exactly 1 time
    21176141 (24.09%) aligned >1 times
98.77% overall alignment rate
[bam_sort_core] merging from 28 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 1:07:57. Running peak memory: 4.632GB.  
  PID: 264831;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 264826;	Command: bowtie2;	Return code: 0;	Memory used: 3.722GB  
  PID: 264832;	Command: samtools;	Return code: 0;	Memory used: 0.891GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam` (271209)
<pre>
</pre>
Command completed. Elapsed time: 0:02:51. Running peak memory: 4.632GB.  
  PID: 271209;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	86810990	PEPPRO	_RES_

> `QC_filtered_reads`	9341224	PEPPRO	_RES_

> `Aligned_reads`	77469766	PEPPRO	_RES_

> `Alignment_rate`	80.01	PEPPRO	_RES_

> `Total_efficiency`	78.01	PEPPRO	_RES_

> `Read_depth`	7.58	PEPPRO	_RES_

### Compress all unmapped read files (06-15 09:07:37) elapsed: 5834.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/prealignments/K562_PRO-seq_20_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/prealignments/K562_PRO-seq_20_human_rDNA_unmap.fq` (274046)
<pre>
</pre>
Command completed. Elapsed time: 0:03:02. Running peak memory: 4.632GB.  
  PID: 274046;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_temp.bam` (274469)
<pre>
</pre>
Command completed. Elapsed time: 0:01:27. Running peak memory: 4.632GB.  
  PID: 274469;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	1829855	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam` (274545)
<pre>
</pre>
Command completed. Elapsed time: 0:01:15. Running peak memory: 4.632GB.  
  PID: 274545;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/chr_sizes.bed` (274607,274608,274609,274610)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 274608;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 274610;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 274607;	Command: samtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 274609;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_noMT.bam` (274612)
<pre>
</pre>
Command completed. Elapsed time: 0:01:42. Running peak memory: 4.632GB.  
  PID: 274612;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam` (274904)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.632GB.  
  PID: 274904;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam` (274906)
<pre>
</pre>
Command completed. Elapsed time: 0:01:15. Running peak memory: 4.632GB.  
  PID: 274906;	Command: samtools;	Return code: 0;	Memory used: 0.013GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 09:20:19) elapsed: 762.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_bamQC.tsv` (275415)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/tmp_K562_PRO-seq_20_sort_mqbelqn0'
Processing with 12 cores...
Discarding 91 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 104 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:01:52. Running peak memory: 4.632GB.  
  PID: 275415;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 2.42GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_bamQC.tsv`

> `NRF`	0.69	PEPPRO	_RES_

> `PBC1`	0.84	PEPPRO	_RES_

> `PBC2`	8.4	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_unmap.bam`  

> `samtools view -b -@ 12 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_unmap.bam` (275537)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 4.632GB.  
  PID: 275537;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_temp.bam`

> `Unmapped_reads`	1079115	PEPPRO	_RES_

### Split BAM by strand (06-15 09:22:43) elapsed: 144.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_plus.bam` (275595)
<pre>
</pre>
Command completed. Elapsed time: 0:06:07. Running peak memory: 4.632GB.  
  PID: 275595;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_minus.bam` (276096)
<pre>
</pre>
Command completed. Elapsed time: 0:05:54. Running peak memory: 4.632GB.  
  PID: 276096;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 09:34:44) elapsed: 721.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (276690)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 276690;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_plus_TssEnrichment.txt` (276695)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 4.632GB.  
  PID: 276695;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.544GB


> `TSS_coding_score`	13.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_minus_TssEnrichment.txt` (276733)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 4.632GB.  
  PID: 276733;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.555GB


> `TSS_non-coding_score`	5.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_minus_TssEnrichment.txt` (276960)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 4.632GB.  
  PID: 276960;	Command: Rscript;	Return code: 0;	Memory used: 0.286GB

> `TSS enrichment`	QC_hg38/K562_PRO-seq_20_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_20_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt` (276986,276987,276988,276989)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 276986;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 276988;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 276987;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 276989;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_keep.txt` (276992)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 276992;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 09:35:16) elapsed: 32.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/hg38_ensembl_tss.bed` (276994,276995)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 4.632GB.  
  PID: 276994;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 276995;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/hg38_ensembl_gene_body.bed` (276998,276999)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 276998;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 276999;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_TSS_density.bed` (277001,277002,277003,277004)
<pre>
</pre>
Command completed. Elapsed time: 0:01:53. Running peak memory: 4.632GB.  
  PID: 277004;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 277001;	Command: bedtools;	Return code: 0;	Memory used: 0.025GB  
  PID: 277003;	Command: sort;	Return code: 0;	Memory used: 0.01GB  
  PID: 277002;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_gene_body_density.bed` (277100,277101,277102)
<pre>
</pre>
Command completed. Elapsed time: 0:02:55. Running peak memory: 4.632GB.  
  PID: 277100;	Command: bedtools;	Return code: 0;	Memory used: 0.224GB  
  PID: 277102;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 277101;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/tmpqfzpdbsd` (277484,277485,277486)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 277484;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 277486;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 277485;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/tmpqfzpdbsd | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0567686) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/tmpqfzpdbsd > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_pause_index.bed` (277492)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 277492;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.75	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_pause_index.bed` (277497)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.632GB.  
  PID: 277497;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Pause index`	QC_hg38/K562_PRO-seq_20_pause_index.pdf	Pause index	QC_hg38/K562_PRO-seq_20_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_pause_index.bed` (277521)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 277521;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 09:40:12) elapsed: 296.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_plus.bam`
77469766 27628563

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_minus.bam`
77469766 26288184

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/signal_hg38/K562_PRO-seq_20_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/hg38_gene_sort.bed` (277654,277655)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.632GB.  
  PID: 277654;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 277655;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/signal_hg38/K562_PRO-seq_20_gene_coverage.bed` (277658)
<pre>
</pre>
Command completed. Elapsed time: 0:02:59. Running peak memory: 4.632GB.  
  PID: 277658;	Command: bedtools;	Return code: 0;	Memory used: 0.222GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/raw/hg38_annotations.bed.gz` (278001)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 278001;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/raw/hg38_annotations.bed` (278002)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 278002;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 09:45:32) elapsed: 320.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/raw/hg38_annotations.bed` (278011)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 4.632GB.  
  PID: 278011;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Enhancer_sort.bed` (278013,278014,278015,278016)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.632GB.  
  PID: 278013;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 278014;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 278016;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB  
  PID: 278015;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Enhancer_plus_coverage.bed` (278019)
<pre>
</pre>
Command completed. Elapsed time: 0:00:56. Running peak memory: 4.632GB.  
  PID: 278019;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Enhancer_minus_coverage.bed` (278067)
<pre>
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 4.632GB.  
  PID: 278067;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter_sort.bed` (278113,278114,278115,278116)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 278113;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 278114;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 278116;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 278115;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Promoter_plus_coverage.bed` (278118)
<pre>
</pre>
Command completed. Elapsed time: 0:01:10. Running peak memory: 4.632GB.  
  PID: 278118;	Command: bedtools;	Return code: 0;	Memory used: 0.115GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Promoter_minus_coverage.bed` (278177)
<pre>
</pre>
Command completed. Elapsed time: 0:01:06. Running peak memory: 4.632GB.  
  PID: 278177;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter_Flanking_Region"` (278233)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 278233;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter_Flanking_Region_sort.bed` (278234,278235,278236,278237)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.632GB.  
  PID: 278234;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 278236;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 278235;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 278237;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Promoter_Flanking_Region_plus_coverage.bed` (278241)
<pre>
</pre>
Command completed. Elapsed time: 0:01:02. Running peak memory: 4.632GB.  
  PID: 278241;	Command: bedtools;	Return code: 0;	Memory used: 0.026GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Promoter_Flanking_Region_minus_coverage.bed` (278534)
<pre>
</pre>
Command completed. Elapsed time: 0:01:03. Running peak memory: 4.632GB.  
  PID: 278534;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/5_UTR"` (278588)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 278588;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/5_UTR_sort.bed` (278589,278590,278591,278592)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.632GB.  
  PID: 278589;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 278590;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 278592;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 278591;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_5_UTR_plus_coverage.bed` (278594)
<pre>
</pre>
Command completed. Elapsed time: 0:00:56. Running peak memory: 4.632GB.  
  PID: 278594;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_5_UTR_minus_coverage.bed` (278643)
<pre>
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 4.632GB.  
  PID: 278643;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/3_UTR"` (278690)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 278690;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/3_UTR_sort.bed` (278691,278692,278693,278694)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.632GB.  
  PID: 278691;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 278692;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 278694;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 278693;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_3_UTR_plus_coverage.bed` (278696)
<pre>
</pre>
Command completed. Elapsed time: 0:01:00. Running peak memory: 4.632GB.  
  PID: 278696;	Command: bedtools;	Return code: 0;	Memory used: 0.026GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_3_UTR_minus_coverage.bed` (278746)
<pre>
</pre>
Command completed. Elapsed time: 0:00:58. Running peak memory: 4.632GB.  
  PID: 278746;	Command: bedtools;	Return code: 0;	Memory used: 0.038GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Exon_sort.bed` (278990,278991,278992,278993)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 4.632GB.  
  PID: 278990;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 278991;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 278993;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 278992;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Exon_plus_coverage.bed` (278997)
<pre>
</pre>
Command completed. Elapsed time: 0:01:03. Running peak memory: 4.632GB.  
  PID: 278997;	Command: bedtools;	Return code: 0;	Memory used: 0.101GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Exon_minus_coverage.bed` (279050)
<pre>
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 4.632GB.  
  PID: 279050;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Intron_sort.bed` (279102,279103,279104,279106)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.632GB.  
  PID: 279102;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 279104;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 279103;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 279106;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Intron_plus_coverage.bed` (279109)
<pre>
</pre>
Command completed. Elapsed time: 0:01:09. Running peak memory: 4.632GB.  
  PID: 279109;	Command: bedtools;	Return code: 0;	Memory used: 0.083GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Intron_minus_coverage.bed` (279168)
<pre>
</pre>
Command completed. Elapsed time: 0:01:10. Running peak memory: 4.632GB.  
  PID: 279168;	Command: bedtools;	Return code: 0;	Memory used: 0.096GB


### Plot cFRiF/FRiF (06-15 10:00:05) elapsed: 872.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_20 -z 3099922541 -n 38543996 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Intron_plus_coverage.bed` (279504)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 4.632GB.  
  PID: 279504;	Command: Rscript;	Return code: 0;	Memory used: 0.439GB

> `cFRiF`	QC_hg38/K562_PRO-seq_20_cFRiF.pdf	cFRiF	QC_hg38/K562_PRO-seq_20_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_20 -z 3099922541 -n 38543996 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_Intron_plus_coverage.bed` (279550)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 4.632GB.  
  PID: 279550;	Command: Rscript;	Return code: 0;	Memory used: 0.439GB

> `FRiF`	QC_hg38/K562_PRO-seq_20_FRiF.pdf	FRiF	QC_hg38/K562_PRO-seq_20_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 10:01:22) elapsed: 78.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/hg38_exons_sort.bed` (279600,279601)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.632GB.  
  PID: 279601;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB  
  PID: 279600;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/hg38_introns_sort.bed` (279608,279609,279610)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.632GB.  
  PID: 279608;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 279610;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 279609;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_exons_coverage.bed` (279616)
<pre>
</pre>
Command completed. Elapsed time: 0:02:03. Running peak memory: 4.632GB.  
  PID: 279616;	Command: bedtools;	Return code: 0;	Memory used: 0.065GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_introns_coverage.bed` (279755)
<pre>
</pre>
Command completed. Elapsed time: 0:02:49. Running peak memory: 4.632GB.  
  PID: 279755;	Command: bedtools;	Return code: 0;	Memory used: 0.11GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/77.469766)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_exons_rpkm.bed` (280091,280092,280093)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.632GB.  
  PID: 280091;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 280093;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 280092;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/77.469766)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_introns_rpkm.bed` (280095,280096,280097)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.632GB.  
  PID: 280095;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 280097;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 280096;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_exon_intron_ratios.bed` (280100,280101,280102)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 280100;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 280102;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 280101;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.39	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_exon_intron_ratios.bed --annotate` (280108)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.632GB.  
  PID: 280108;	Command: Rscript;	Return code: 0;	Memory used: 0.285GB

> `mRNA contamination`	QC_hg38/K562_PRO-seq_20_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_PRO-seq_20_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/QC_hg38/K562_PRO-seq_20_exon_intron_ratios.bed` (280129)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.632GB.  
  PID: 280129;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 10:06:32) elapsed: 310.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/signal_hg38/K562_PRO-seq_20_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/signal_hg38/K562_PRO-seq_20_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_plus.bam` (280138)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 4.632GB.  
  PID: 280138;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/signal_hg38/K562_PRO-seq_20_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/signal_hg38/K562_PRO-seq_20_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 77469766.0` (280170)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_plus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_20_plus_cuttrace_wlhp8s6l'
Processing with 4 cores...
stdin is empty of data
Discarding 101 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr14_KI270726v1_random', 'chr22_KI270732v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 94 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 94 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/signal_hg38/K562_PRO-seq_20_plus_exact_body_0-mer.bw'
Merging 94 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/signal_hg38/K562_PRO-seq_20_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:07. Running peak memory: 4.632GB.  
  PID: 280170;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.465GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/signal_hg38/K562_PRO-seq_20_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/signal_hg38/K562_PRO-seq_20_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_minus.bam` (281940)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 4.632GB.  
  PID: 281940;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/signal_hg38/K562_PRO-seq_20_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/signal_hg38/K562_PRO-seq_20_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 77469766.0` (281970)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/aligned_hg38/K562_PRO-seq_20_minus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_20_minus_cuttrace_uzz4vhui'
Processing with 4 cores...
stdin is empty of data
stdin is empty of data
Discarding 105 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 90 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 90 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/signal_hg38/K562_PRO-seq_20_minus_exact_body_0-mer.bw'
Merging 90 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_20/signal_hg38/K562_PRO-seq_20_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:03. Running peak memory: 4.632GB.  
  PID: 281970;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.19GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  3:28:03
*  Total elapsed time (all runs):  4:40:15
*         Peak memory (this run):  4.6316 GB
*        Pipeline completed time: 2020-06-15 10:25:55
