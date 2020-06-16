### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_100 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_100pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba26-33c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/
*  Pipeline started at:   (06-15 07:17:12) elapsed: 1.0 _TIME_

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
*              `cores`:  `12`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_100pctRNA.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_100`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_100pctRNA.fastq.gz

> `File_mb`	5091.82	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:13) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/raw/K562_RNA-seq_100.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_100pctRNA.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/raw/K562_RNA-seq_100.fastq.gz` (52656)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 52656;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/raw/K562_RNA-seq_100.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastq/K562_RNA-seq_100_R1.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/raw/K562_RNA-seq_100.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastq/K562_RNA-seq_100_R1.fastq` (52657)
<pre>
</pre>
Command completed. Elapsed time: 0:02:35. Running peak memory: 0.002GB.  
  PID: 52657;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	70000000	PEPPRO	_RES_

> `Fastq_reads`	70000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/raw/K562_RNA-seq_100.fastq.gz']

### FASTQ processing:  (06-15 07:21:10) elapsed: 237.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastq/K562_RNA-seq_100_R1_processed.fastq`  

> `(cutadapt -j 12 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastq/K562_RNA-seq_100_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastq/K562_RNA-seq_100_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastq/K562_RNA-seq_100_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/cutadapt/K562_RNA-seq_100_R1_cutadapt.txt` (53122)
<pre>
</pre>
Command completed. Elapsed time: 0:01:53. Running peak memory: 5.727GB.  
  PID: 53122;	Command: cutadapt;	Return code: 0;	Memory used: 5.727GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastq/K562_RNA-seq_100_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastq/K562_RNA-seq_100_R1_processed.fastq` (53248,53249)
<pre>
</pre>
Command completed. Elapsed time: 0:01:58. Running peak memory: 5.727GB.  
  PID: 53248;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 53249;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/cutadapt/K562_RNA-seq_100_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	24858626.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/cutadapt/K562_RNA-seq_100_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/cutadapt/K562_RNA-seq_100_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastq/K562_RNA-seq_100_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	0.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	0.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastqc/K562_RNA-seq_100_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (53615)
### Calculate the number of trimmed reads
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 53615;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	70000000	PEPPRO	_RES_

> `Trim_loss_rate`	0.0	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastq/K562_RNA-seq_100_R1_processed.fastq` (53636)
<pre>
Started analysis of K562_RNA-seq_100_R1_processed.fastq
Approx 5% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 10% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 15% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 20% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 25% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 30% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 35% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 40% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 45% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 50% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 55% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 60% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 65% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 70% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 75% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 80% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 85% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 90% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 95% complete for K562_RNA-seq_100_R1_processed.fastq
Approx 100% complete for K562_RNA-seq_100_R1_processed.fastq
Analysis complete for K562_RNA-seq_100_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:04:35. Running peak memory: 5.727GB.  
  PID: 53636;	Command: fastqc;	Return code: 0;	Memory used: 0.201GB

> `FastQC report r1`	fastqc/K562_RNA-seq_100_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastq/processed_R1.flag` (54166)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 54166;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-15 07:30:55) elapsed: 585.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/cutadapt/K562_RNA-seq_100_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/cutadapt/K562_RNA-seq_100_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/cutadapt` (54167)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.727GB.  
  PID: 54167;	Command: Rscript;	Return code: 0;	Memory used: 0.3GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_100_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_100_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/cutadapt/K562_RNA-seq_100_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	0	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:31:00) elapsed: 5.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/cutadapt/K562_RNA-seq_100_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/cutadapt/K562_RNA-seq_100_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/cutadapt/K562_RNA-seq_100_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/cutadapt/K562_RNA-seq_100_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/cutadapt/K562_RNA-seq_100_R1_cutadapt.txt | awk 'END {print $1}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/cutadapt/K562_RNA-seq_100_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 10 && $1 >= 16){degradedSum += $2}; ($1 >= 6 && $1 <= 0){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.0	PEPPRO	_RES_

### Prealignments (06-15 07:31:00) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:31:00) elapsed: 0.0 _TIME_


> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_100 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/fastq/K562_RNA-seq_100_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/prealignments/K562_RNA-seq_100_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
70000000 reads; of these:
  70000000 (100.00%) were unpaired; of these:
    68375974 (97.68%) aligned 0 times
    1624026 (2.32%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
2.32% overall alignment rate

> `Aligned_reads_human_rDNA`	1624026.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	2.32	PEPPRO	_RES_

### Map to genome (06-15 07:39:19) elapsed: 500.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam`  

> `bowtie2 -p 12 --very-sensitive --rg-id K562_RNA-seq_100 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/prealignments/K562_RNA-seq_100_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/tmpd6hrvyrl -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_temp.bam` (54845,54852,54853)
<pre>
68375974 reads; of these:
  68375974 (100.00%) were unpaired; of these:
    11525030 (16.86%) aligned 0 times
    33591968 (49.13%) aligned exactly 1 time
    23258976 (34.02%) aligned >1 times
83.14% overall alignment rate
[bam_sort_core] merging from 22 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:36:44. Running peak memory: 5.727GB.  
  PID: 54845;	Command: bowtie2;	Return code: 0;	Memory used: 3.692GB  
  PID: 54852;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 54853;	Command: samtools;	Return code: 0;	Memory used: 0.89GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam` (58535)
<pre>
</pre>
Command completed. Elapsed time: 0:03:54. Running peak memory: 5.727GB.  
  PID: 58535;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	56850944	PEPPRO	_RES_

> `QC_filtered_reads`	12689564	PEPPRO	_RES_

> `Aligned_reads`	44161380	PEPPRO	_RES_

> `Alignment_rate`	63.09	PEPPRO	_RES_

> `Total_efficiency`	63.09	PEPPRO	_RES_

> `Read_depth`	11.83	PEPPRO	_RES_

### Compress all unmapped read files (06-15 08:32:47) elapsed: 3208.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/prealignments/K562_RNA-seq_100_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/prealignments/K562_RNA-seq_100_human_rDNA_unmap.fq` (60113)
<pre>
</pre>
Command completed. Elapsed time: 0:02:33. Running peak memory: 5.727GB.  
  PID: 60113;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_temp.bam` (60471)
<pre>
</pre>
Command completed. Elapsed time: 0:01:14. Running peak memory: 5.727GB.  
  PID: 60471;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	4192362	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam` (60538)
<pre>
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 5.727GB.  
  PID: 60538;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/chr_sizes.bed` (60579,60580,60581,60582)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 60579;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 60581;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 60580;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 60582;	Command: grep;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_noMT.bam` (60584)
<pre>
</pre>
Command completed. Elapsed time: 0:00:53. Running peak memory: 5.727GB.  
  PID: 60584;	Command: samtools;	Return code: 0;	Memory used: 0.02GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam` (60644)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 60644;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam` (60645)
<pre>
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 5.727GB.  
  PID: 60645;	Command: samtools;	Return code: 0;	Memory used: 0.014GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	76	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 08:41:13) elapsed: 505.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_bamQC.tsv` (61041)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/tmp_K562_RNA-seq_100_sort_iqht1p80'
Processing with 12 cores...
Discarding 104 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270723v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 91 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:48. Running peak memory: 5.727GB.  
  PID: 61041;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.987GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_bamQC.tsv`

> `NRF`	0.47	PEPPRO	_RES_

> `PBC1`	0.79	PEPPRO	_RES_

> `PBC2`	8.82	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_unmap.bam`  

> `samtools view -b -@ 12 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_unmap.bam` (61111)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 5.727GB.  
  PID: 61111;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_temp.bam`

> `Unmapped_reads`	11525030	PEPPRO	_RES_

### Split BAM by strand (06-15 08:42:36) elapsed: 83.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_plus.bam` (61169)
<pre>
</pre>
Command completed. Elapsed time: 0:03:43. Running peak memory: 5.727GB.  
  PID: 61169;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_minus.bam` (61550)
<pre>
</pre>
Command completed. Elapsed time: 0:03:22. Running peak memory: 5.727GB.  
  PID: 61550;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 08:49:41) elapsed: 424.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (61718)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 61718;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_plus_TssEnrichment.txt` (61719)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 5.727GB.  
  PID: 61719;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.515GB


> `TSS_coding_score`	19.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_minus_TssEnrichment.txt` (61754)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 5.727GB.  
  PID: 61754;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.59GB


> `TSS_non-coding_score`	2.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_minus_TssEnrichment.txt` (61787)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 5.727GB.  
  PID: 61787;	Command: Rscript;	Return code: 0;	Memory used: 0.308GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_100_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_100_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt` (62048,62049,62050,62051)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 62048;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 62050;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 62049;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 62051;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_keep.txt` (62054)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 62054;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-15 08:50:05) elapsed: 24.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/hg38_ensembl_tss.bed` (62056,62057)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 5.727GB.  
  PID: 62056;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 62057;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/hg38_ensembl_gene_body.bed` (62060,62061)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 62060;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 62061;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_TSS_density.bed` (62063,62064,62065,62066)
<pre>
</pre>
Command completed. Elapsed time: 0:01:08. Running peak memory: 5.727GB.  
  PID: 62063;	Command: bedtools;	Return code: 0;	Memory used: 0.072GB  
  PID: 62065;	Command: sort;	Return code: 0;	Memory used: 0.01GB  
  PID: 62064;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 62066;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_gene_body_density.bed` (62127,62128,62129)
<pre>
</pre>
Command completed. Elapsed time: 0:01:46. Running peak memory: 5.727GB.  
  PID: 62127;	Command: bedtools;	Return code: 0;	Memory used: 0.313GB  
  PID: 62129;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 62128;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/tmpv201glls` (62218,62219,62220)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 62218;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 62220;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 62219;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/tmpv201glls | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0486842) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/tmpv201glls > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_pause_index.bed` (62226)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 62226;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	7.79	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_pause_index.bed` (62232)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.727GB.  
  PID: 62232;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Pause index`	QC_hg38/K562_RNA-seq_100_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_100_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_pause_index.bed` (62253)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 62253;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 08:53:07) elapsed: 183.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_plus.bam`
44161380 17667752

> `Plus_FRiP`	0.4	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_minus.bam`
44161380 17985866

> `Minus_FRiP`	0.41	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/signal_hg38/K562_RNA-seq_100_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/hg38_gene_sort.bed` (62331,62332)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.727GB.  
  PID: 62331;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 62332;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/signal_hg38/K562_RNA-seq_100_gene_coverage.bed` (62335)
<pre>
</pre>
Command completed. Elapsed time: 0:01:41. Running peak memory: 5.727GB.  
  PID: 62335;	Command: bedtools;	Return code: 0;	Memory used: 0.269GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/raw/hg38_annotations.bed.gz` (62629)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 62629;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/raw/hg38_annotations.bed` (62630)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 62630;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 08:56:05) elapsed: 177.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/raw/hg38_annotations.bed` (62639)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.727GB.  
  PID: 62639;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Enhancer_sort.bed` (62641,62642,62643,62644)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.727GB.  
  PID: 62641;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 62642;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 62644;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 62643;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Enhancer_plus_coverage.bed` (62647)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 5.727GB.  
  PID: 62647;	Command: bedtools;	Return code: 0;	Memory used: 0.038GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Enhancer_minus_coverage.bed` (62676)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 5.727GB.  
  PID: 62676;	Command: bedtools;	Return code: 0;	Memory used: 0.049GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter_sort.bed` (62706,62707,62708,62709)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 62706;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 62707;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 62709;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 62708;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Promoter_plus_coverage.bed` (62711)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 5.727GB.  
  PID: 62711;	Command: bedtools;	Return code: 0;	Memory used: 0.195GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Promoter_minus_coverage.bed` (62743)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 5.727GB.  
  PID: 62743;	Command: bedtools;	Return code: 0;	Memory used: 0.09GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter_Flanking_Region"` (62776)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 62776;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter_Flanking_Region_sort.bed` (62777,62778,62779,62780)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.727GB.  
  PID: 62777;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 62779;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 62778;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 62780;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Promoter_Flanking_Region_plus_coverage.bed` (62783)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 5.727GB.  
  PID: 62783;	Command: bedtools;	Return code: 0;	Memory used: 0.026GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Promoter_Flanking_Region_minus_coverage.bed` (62812)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 5.727GB.  
  PID: 62812;	Command: bedtools;	Return code: 0;	Memory used: 0.069GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/5_UTR"` (62841)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 62841;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/5_UTR_sort.bed` (62842,62843,62844,62845)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.727GB.  
  PID: 62842;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 62843;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 62845;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 62844;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_5_UTR_plus_coverage.bed` (62848)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 5.727GB.  
  PID: 62848;	Command: bedtools;	Return code: 0;	Memory used: 0.031GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_5_UTR_minus_coverage.bed` (63134)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 5.727GB.  
  PID: 63134;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/3_UTR"` (63196)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 63196;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/3_UTR_sort.bed` (63197,63198,63199,63200)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.727GB.  
  PID: 63197;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 63198;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 63200;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 63199;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_3_UTR_plus_coverage.bed` (63203)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 5.727GB.  
  PID: 63203;	Command: bedtools;	Return code: 0;	Memory used: 0.06GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_3_UTR_minus_coverage.bed` (63250)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 5.727GB.  
  PID: 63250;	Command: bedtools;	Return code: 0;	Memory used: 0.076GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Exon_sort.bed` (63284,63285,63286,63287)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 5.727GB.  
  PID: 63284;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 63285;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 63287;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 63286;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Exon_plus_coverage.bed` (63291)
<pre>
</pre>
Command completed. Elapsed time: 0:00:48. Running peak memory: 5.727GB.  
  PID: 63291;	Command: bedtools;	Return code: 0;	Memory used: 0.248GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Exon_minus_coverage.bed` (63332)
<pre>
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 5.727GB.  
  PID: 63332;	Command: bedtools;	Return code: 0;	Memory used: 0.132GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Intron_sort.bed` (63372,63373,63374,63375)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 5.727GB.  
  PID: 63372;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 63374;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 63373;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 63375;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Intron_plus_coverage.bed` (63378)
<pre>
</pre>
Command completed. Elapsed time: 0:00:43. Running peak memory: 5.727GB.  
  PID: 63378;	Command: bedtools;	Return code: 0;	Memory used: 0.18GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Intron_minus_coverage.bed` (63414)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 5.727GB.  
  PID: 63414;	Command: bedtools;	Return code: 0;	Memory used: 0.089GB


### Plot cFRiF/FRiF (06-15 09:05:01) elapsed: 536.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_RNA-seq_100 -z 3099922541 -n 20294129 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Intron_plus_coverage.bed` (63657)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 5.727GB.  
  PID: 63657;	Command: Rscript;	Return code: 0;	Memory used: 0.46GB

> `cFRiF`	QC_hg38/K562_RNA-seq_100_cFRiF.pdf	cFRiF	QC_hg38/K562_RNA-seq_100_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_RNA-seq_100 -z 3099922541 -n 20294129 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_Intron_plus_coverage.bed` (63701)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 5.727GB.  
  PID: 63701;	Command: Rscript;	Return code: 0;	Memory used: 0.442GB

> `FRiF`	QC_hg38/K562_RNA-seq_100_FRiF.pdf	FRiF	QC_hg38/K562_RNA-seq_100_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 09:06:10) elapsed: 69.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/hg38_exons_sort.bed` (63738,63739)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.727GB.  
  PID: 63739;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB  
  PID: 63738;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/hg38_introns_sort.bed` (63745,63746,63747)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.727GB.  
  PID: 63745;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 63747;	Command: bedtools;	Return code: 0;	Memory used: 0.003GB  
  PID: 63746;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_exons_coverage.bed` (63754)
<pre>
</pre>
Command completed. Elapsed time: 0:01:30. Running peak memory: 5.727GB.  
  PID: 63754;	Command: bedtools;	Return code: 0;	Memory used: 0.22GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_introns_coverage.bed` (63829)
<pre>
</pre>
Command completed. Elapsed time: 0:01:28. Running peak memory: 5.727GB.  
  PID: 63829;	Command: bedtools;	Return code: 0;	Memory used: 0.206GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/44.16138)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_exons_rpkm.bed` (63902,63903,63904)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.727GB.  
  PID: 63902;	Command: awk;	Return code: 0;	Memory used: 0.006GB  
  PID: 63904;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 63903;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/44.16138)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_introns_rpkm.bed` (63907,63908,63909)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.727GB.  
  PID: 63907;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 63909;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 63908;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_exon_intron_ratios.bed` (63911,63912,63913)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 63911;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 63913;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 63912;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	12.77	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_exon_intron_ratios.bed --annotate` (63919)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.727GB.  
  PID: 63919;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_100_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_100_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/QC_hg38/K562_RNA-seq_100_exon_intron_ratios.bed` (63940)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.727GB.  
  PID: 63940;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 09:09:25) elapsed: 195.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/signal_hg38/K562_RNA-seq_100_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/signal_hg38/K562_RNA-seq_100_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_plus.bam` (63948)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 5.727GB.  
  PID: 63948;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/signal_hg38/K562_RNA-seq_100_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/signal_hg38/K562_RNA-seq_100_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 44161380.0` (63968)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_100_plus_cuttrace_i5cpm31b'
Processing with 4 cores...
Discarding 116 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1']
Keeping 79 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270726v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 79 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/signal_hg38/K562_RNA-seq_100_plus_exact_body_0-mer.bw'
Merging 79 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/signal_hg38/K562_RNA-seq_100_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:08. Running peak memory: 5.727GB.  
  PID: 63968;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.986GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/signal_hg38/K562_RNA-seq_100_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/signal_hg38/K562_RNA-seq_100_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_minus.bam` (65502)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 5.727GB.  
  PID: 65502;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/signal_hg38/K562_RNA-seq_100_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/signal_hg38/K562_RNA-seq_100_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 44161380.0` (65522)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/aligned_hg38/K562_RNA-seq_100_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_100_minus_cuttrace_za4nvju3'
Processing with 4 cores...
Discarding 118 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 77 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270732v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270583v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 77 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/signal_hg38/K562_RNA-seq_100_minus_exact_body_0-mer.bw'
Merging 77 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_100/signal_hg38/K562_RNA-seq_100_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:08. Running peak memory: 5.727GB.  
  PID: 65522;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.701GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  2:07:14
*  Total elapsed time (all runs):  2:42:53
*         Peak memory (this run):  5.7268 GB
*        Pipeline completed time: 2020-06-15 09:24:25
