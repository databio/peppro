### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_0 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_0pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aw29-26a
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/
*  Pipeline started at:   (06-15 07:17:12) elapsed: 5.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_0pctRNA.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_0`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_0pctRNA.fastq.gz

> `File_mb`	5665.7	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:13) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/raw/K562_RNA-seq_0.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_0pctRNA.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/raw/K562_RNA-seq_0.fastq.gz` (436786)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 436786;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/raw/K562_RNA-seq_0.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastq/K562_RNA-seq_0_R1.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/raw/K562_RNA-seq_0.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastq/K562_RNA-seq_0_R1.fastq` (436787)
<pre>
</pre>
Command completed. Elapsed time: 0:03:09. Running peak memory: 0.002GB.  
  PID: 436787;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	70000000	PEPPRO	_RES_

> `Fastq_reads`	70000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/raw/K562_RNA-seq_0.fastq.gz']

### FASTQ processing:  (06-15 07:22:04) elapsed: 292.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastq/K562_RNA-seq_0_R1_processed.fastq`  

> `(cutadapt -j 12 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastq/K562_RNA-seq_0_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastq/K562_RNA-seq_0_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastq/K562_RNA-seq_0_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/cutadapt/K562_RNA-seq_0_R1_cutadapt.txt` (437302)
<pre>
</pre>
Command completed. Elapsed time: 0:02:15. Running peak memory: 4.588GB.  
  PID: 437302;	Command: cutadapt;	Return code: 0;	Memory used: 4.588GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastq/K562_RNA-seq_0_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastq/K562_RNA-seq_0_R1_processed.fastq` (437447,437448)
<pre>
</pre>
Command completed. Elapsed time: 0:02:12. Running peak memory: 4.588GB.  
  PID: 437448;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB  
  PID: 437447;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/cutadapt/K562_RNA-seq_0_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	59632858.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/cutadapt/K562_RNA-seq_0_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/cutadapt/K562_RNA-seq_0_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastq/K562_RNA-seq_0_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	1752284.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	2.5033	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastqc/K562_RNA-seq_0_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (437831)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 437831;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	68247716	PEPPRO	_RES_

> `Trim_loss_rate`	2.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastq/K562_RNA-seq_0_R1_processed.fastq` (437855)
<pre>
Started analysis of K562_RNA-seq_0_R1_processed.fastq
Approx 5% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 10% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 15% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 20% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 25% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 30% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 35% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 40% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 45% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 50% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 55% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 60% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 65% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 70% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 75% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 80% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 85% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 90% complete for K562_RNA-seq_0_R1_processed.fastq
Approx 95% complete for K562_RNA-seq_0_R1_processed.fastq
Analysis complete for K562_RNA-seq_0_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:04:44. Running peak memory: 4.588GB.  
  PID: 437855;	Command: fastqc;	Return code: 0;	Memory used: 0.21GB

> `FastQC report r1`	fastqc/K562_RNA-seq_0_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastq/processed_R1.flag` (438396)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 438396;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-15 07:32:44) elapsed: 640.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/cutadapt/K562_RNA-seq_0_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/cutadapt/K562_RNA-seq_0_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/cutadapt` (438397)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 4.588GB.  
  PID: 438397;	Command: Rscript;	Return code: 0;	Memory used: 0.116GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_0_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_0_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/cutadapt/K562_RNA-seq_0_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:32:53) elapsed: 9.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/cutadapt/K562_RNA-seq_0_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/cutadapt/K562_RNA-seq_0_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/cutadapt/K562_RNA-seq_0_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/cutadapt/K562_RNA-seq_0_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/cutadapt/K562_RNA-seq_0_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2311	PEPPRO	_RES_

### Prealignments (06-15 07:32:53) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:32:53) elapsed: 0.0 _TIME_


> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_0 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/fastq/K562_RNA-seq_0_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/prealignments/K562_RNA-seq_0_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
68247716 reads; of these:
  68247716 (100.00%) were unpaired; of these:
    61949206 (90.77%) aligned 0 times
    6298510 (9.23%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.23% overall alignment rate

> `Aligned_reads_human_rDNA`	6298510.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.23	PEPPRO	_RES_

### Map to genome (06-15 07:39:59) elapsed: 426.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam`  

> `bowtie2 -p 12 --very-sensitive --rg-id K562_RNA-seq_0 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/prealignments/K562_RNA-seq_0_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/tmpvvqfo_w7 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_temp.bam` (439014,439020,439021)
<pre>
61949206 reads; of these:
  61949206 (100.00%) were unpaired; of these:
    760871 (1.23%) aligned 0 times
    46260822 (74.68%) aligned exactly 1 time
    14927513 (24.10%) aligned >1 times
98.77% overall alignment rate
[bam_sort_core] merging from 19 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:47:12. Running peak memory: 4.588GB.  
  PID: 439020;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 439014;	Command: bowtie2;	Return code: 0;	Memory used: 3.721GB  
  PID: 439021;	Command: samtools;	Return code: 0;	Memory used: 0.889GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam` (443653)
<pre>
</pre>
Command completed. Elapsed time: 0:01:48. Running peak memory: 4.588GB.  
  PID: 443653;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	61188335	PEPPRO	_RES_

> `QC_filtered_reads`	6585652	PEPPRO	_RES_

> `Aligned_reads`	54602683	PEPPRO	_RES_

> `Alignment_rate`	80.01	PEPPRO	_RES_

> `Total_efficiency`	78.0	PEPPRO	_RES_

> `Read_depth`	5.84	PEPPRO	_RES_

### Compress all unmapped read files (06-15 08:51:05) elapsed: 4266.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/prealignments/K562_RNA-seq_0_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/prealignments/K562_RNA-seq_0_human_rDNA_unmap.fq` (446024)
<pre>
</pre>
Command completed. Elapsed time: 0:02:09. Running peak memory: 4.588GB.  
  PID: 446024;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_temp.bam` (446160)
<pre>
</pre>
Command completed. Elapsed time: 0:01:00. Running peak memory: 4.588GB.  
  PID: 446160;	Command: samtools;	Return code: 0;	Memory used: 0.017GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	1289952	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam` (446213)
<pre>
</pre>
Command completed. Elapsed time: 0:00:53. Running peak memory: 4.588GB.  
  PID: 446213;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/chr_sizes.bed` (446447,446448,446449,446450)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 446449;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 446447;	Command: samtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 446450;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 446448;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_noMT.bam` (446452)
<pre>
</pre>
Command completed. Elapsed time: 0:01:06. Running peak memory: 4.588GB.  
  PID: 446452;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam` (446526)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 446526;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam` (446527)
<pre>
</pre>
Command completed. Elapsed time: 0:00:52. Running peak memory: 4.588GB.  
  PID: 446527;	Command: samtools;	Return code: 0;	Memory used: 0.016GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 08:59:56) elapsed: 532.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_bamQC.tsv` (446717)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/tmp_K562_RNA-seq_0_sort_fwh_gzl4'
Processing with 12 cores...
Discarding 96 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270725v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1']
Keeping 99 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:01:17. Running peak memory: 4.588GB.  
  PID: 446717;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.766GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_bamQC.tsv`

> `NRF`	0.76	PEPPRO	_RES_

> `PBC1`	0.88	PEPPRO	_RES_

> `PBC2`	10.64	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_unmap.bam`  

> `samtools view -b -@ 12 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_unmap.bam` (447078)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 4.588GB.  
  PID: 447078;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_temp.bam`

> `Unmapped_reads`	760871	PEPPRO	_RES_

### Split BAM by strand (06-15 09:01:35) elapsed: 99.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_plus.bam` (447126)
<pre>
</pre>
Command completed. Elapsed time: 0:04:15. Running peak memory: 4.588GB.  
  PID: 447126;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_minus.bam` (447568)
<pre>
</pre>
Command completed. Elapsed time: 0:04:03. Running peak memory: 4.588GB.  
  PID: 447568;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 09:09:53) elapsed: 498.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (447769)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 447769;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_plus_TssEnrichment.txt` (447771)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 4.588GB.  
  PID: 447771;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.536GB


> `TSS_coding_score`	13.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_minus_TssEnrichment.txt` (448048)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 4.588GB.  
  PID: 448048;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.546GB


> `TSS_non-coding_score`	5.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_minus_TssEnrichment.txt` (448085)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 4.588GB.  
  PID: 448085;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_0_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_0_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt` (448107,448108,448109,448110)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 448107;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 448109;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 448108;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 448110;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_keep.txt` (448112)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 448112;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-15 09:10:19) elapsed: 26.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/hg38_ensembl_tss.bed` (448114,448115)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 4.588GB.  
  PID: 448114;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 448115;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/hg38_ensembl_gene_body.bed` (448119,448120)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 448119;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 448120;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_TSS_density.bed` (448122,448123,448124,448125)
<pre>
</pre>
Command completed. Elapsed time: 0:01:20. Running peak memory: 4.588GB.  
  PID: 448123;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 448125;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 448122;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB  
  PID: 448124;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_gene_body_density.bed` (448192,448193,448194)
<pre>
</pre>
Command completed. Elapsed time: 0:02:04. Running peak memory: 4.588GB.  
  PID: 448193;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 448192;	Command: bedtools;	Return code: 0;	Memory used: 0.164GB  
  PID: 448194;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/tmpdworpk23` (448299,448300,448301)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 448299;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 448301;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 448300;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/tmpdworpk23 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.042526) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/tmpdworpk23 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_pause_index.bed` (448307)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 448307;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.69	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_pause_index.bed` (448312)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.588GB.  
  PID: 448312;	Command: Rscript;	Return code: 0;	Memory used: 0.285GB

> `Pause index`	QC_hg38/K562_RNA-seq_0_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_0_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_pause_index.bed` (448333)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 448333;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 09:13:52) elapsed: 213.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_plus.bam`
54602683 19470784

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_minus.bam`
54602683 18528254

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/signal_hg38/K562_RNA-seq_0_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/hg38_gene_sort.bed` (448624,448625)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.588GB.  
  PID: 448624;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 448625;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/signal_hg38/K562_RNA-seq_0_gene_coverage.bed` (448628)
<pre>
</pre>
Command completed. Elapsed time: 0:02:03. Running peak memory: 4.588GB.  
  PID: 448628;	Command: bedtools;	Return code: 0;	Memory used: 0.163GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/raw/hg38_annotations.bed.gz` (448730)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 448730;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/raw/hg38_annotations.bed` (448731)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 448731;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 09:17:31) elapsed: 220.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/raw/hg38_annotations.bed` (448741)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.588GB.  
  PID: 448741;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Enhancer_sort.bed` (448743,448744,448745,448746)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.588GB.  
  PID: 448743;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 448744;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 448746;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB  
  PID: 448745;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Enhancer_plus_coverage.bed` (448748)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 4.588GB.  
  PID: 448748;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Enhancer_minus_coverage.bed` (448783)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 4.588GB.  
  PID: 448783;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter_sort.bed` (448814,448815,448816,448817)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 448814;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 448815;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 448817;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 448816;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Promoter_plus_coverage.bed` (448820)
<pre>
</pre>
Command completed. Elapsed time: 0:00:45. Running peak memory: 4.588GB.  
  PID: 448820;	Command: bedtools;	Return code: 0;	Memory used: 0.089GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Promoter_minus_coverage.bed` (448858)
<pre>
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 4.588GB.  
  PID: 448858;	Command: bedtools;	Return code: 0;	Memory used: 0.072GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter_Flanking_Region"` (449137)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 449137;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter_Flanking_Region_sort.bed` (449139,449140,449141,449142)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.588GB.  
  PID: 449139;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 449141;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 449140;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 449142;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Promoter_Flanking_Region_plus_coverage.bed` (449144)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 4.588GB.  
  PID: 449144;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Promoter_Flanking_Region_minus_coverage.bed` (449180)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 4.588GB.  
  PID: 449180;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/5_UTR"` (449213)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 449213;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/5_UTR_sort.bed` (449214,449215,449216,449217)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.588GB.  
  PID: 449214;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 449215;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 449217;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 449216;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_5_UTR_plus_coverage.bed` (449219)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 4.588GB.  
  PID: 449219;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_5_UTR_minus_coverage.bed` (449255)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 4.588GB.  
  PID: 449255;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/3_UTR"` (449286)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 449286;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/3_UTR_sort.bed` (449287,449288,449289,449290)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.588GB.  
  PID: 449287;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 449288;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 449290;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 449289;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_3_UTR_plus_coverage.bed` (449293)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 4.588GB.  
  PID: 449293;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_3_UTR_minus_coverage.bed` (449328)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 4.588GB.  
  PID: 449328;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Exon_sort.bed` (449363,449364,449365,449366)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 4.588GB.  
  PID: 449363;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 449364;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 449366;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 449365;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Exon_plus_coverage.bed` (449370)
<pre>
</pre>
Command completed. Elapsed time: 0:00:47. Running peak memory: 4.588GB.  
  PID: 449370;	Command: bedtools;	Return code: 0;	Memory used: 0.075GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Exon_minus_coverage.bed` (449603)
<pre>
</pre>
Command completed. Elapsed time: 0:00:44. Running peak memory: 4.588GB.  
  PID: 449603;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Intron_sort.bed` (449644,449645,449646,449647)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.588GB.  
  PID: 449644;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 449646;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 449645;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 449647;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Intron_plus_coverage.bed` (449650)
<pre>
</pre>
Command completed. Elapsed time: 0:00:53. Running peak memory: 4.588GB.  
  PID: 449650;	Command: bedtools;	Return code: 0;	Memory used: 0.076GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Intron_minus_coverage.bed` (449694)
<pre>
</pre>
Command completed. Elapsed time: 0:00:50. Running peak memory: 4.588GB.  
  PID: 449694;	Command: bedtools;	Return code: 0;	Memory used: 0.074GB


### Plot cFRiF/FRiF (06-15 09:27:39) elapsed: 608.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_RNA-seq_0 -z 3099922541 -n 27166534 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Intron_plus_coverage.bed` (449755)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:43. Running peak memory: 4.588GB.  
  PID: 449755;	Command: Rscript;	Return code: 0;	Memory used: 0.446GB

> `cFRiF`	QC_hg38/K562_RNA-seq_0_cFRiF.pdf	cFRiF	QC_hg38/K562_RNA-seq_0_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_RNA-seq_0 -z 3099922541 -n 27166534 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_Intron_plus_coverage.bed` (449802)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 4.588GB.  
  PID: 449802;	Command: Rscript;	Return code: 0;	Memory used: 0.446GB

> `FRiF`	QC_hg38/K562_RNA-seq_0_FRiF.pdf	FRiF	QC_hg38/K562_RNA-seq_0_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 09:28:57) elapsed: 78.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/hg38_exons_sort.bed` (449839,449840)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.588GB.  
  PID: 449840;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB  
  PID: 449839;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/hg38_introns_sort.bed` (449847,449848,449849)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.588GB.  
  PID: 449847;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 449849;	Command: bedtools;	Return code: 0;	Memory used: 0.003GB  
  PID: 449848;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_exons_coverage.bed` (449856)
<pre>
</pre>
Command completed. Elapsed time: 0:01:24. Running peak memory: 4.588GB.  
  PID: 449856;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_introns_coverage.bed` (450182)
<pre>
</pre>
Command completed. Elapsed time: 0:01:53. Running peak memory: 4.588GB.  
  PID: 450182;	Command: bedtools;	Return code: 0;	Memory used: 0.081GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/54.602683)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_exons_rpkm.bed` (450314,450315,450316)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.588GB.  
  PID: 450314;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 450316;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 450315;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/54.602683)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_introns_rpkm.bed` (450318,450319,450320)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.588GB.  
  PID: 450318;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 450320;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 450319;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_exon_intron_ratios.bed` (450323,450324,450325)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 450323;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 450325;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 450324;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.4	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_exon_intron_ratios.bed --annotate` (450331)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.588GB.  
  PID: 450331;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_0_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_0_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/QC_hg38/K562_RNA-seq_0_exon_intron_ratios.bed` (450352)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.588GB.  
  PID: 450352;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 09:32:32) elapsed: 215.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/signal_hg38/K562_RNA-seq_0_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/signal_hg38/K562_RNA-seq_0_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_plus.bam` (450361)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 4.588GB.  
  PID: 450361;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/signal_hg38/K562_RNA-seq_0_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/signal_hg38/K562_RNA-seq_0_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 54602683.0` (450383)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_0_plus_cuttrace_753lphbw'
Processing with 4 cores...
stdin is empty of data
stdin is empty of data
Discarding 109 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270732v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 86 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/signal_hg38/K562_RNA-seq_0_plus_exact_body_0-mer.bw'
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/signal_hg38/K562_RNA-seq_0_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:41. Running peak memory: 4.588GB.  
  PID: 450383;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.902GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/signal_hg38/K562_RNA-seq_0_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/signal_hg38/K562_RNA-seq_0_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_minus.bam` (452055)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 4.588GB.  
  PID: 452055;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/signal_hg38/K562_RNA-seq_0_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/signal_hg38/K562_RNA-seq_0_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 54602683.0` (452078)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/aligned_hg38/K562_RNA-seq_0_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_0_minus_cuttrace_s_3nirwy'
Processing with 4 cores...
stdin is empty of data
stdin is empty of data
Discarding 107 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 88 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 88 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/signal_hg38/K562_RNA-seq_0_minus_exact_body_0-mer.bw'
Merging 88 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_0/signal_hg38/K562_RNA-seq_0_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:30. Running peak memory: 4.588GB.  
  PID: 452078;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.352GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  2:33:27
*  Total elapsed time (all runs):  3:23:42
*         Peak memory (this run):  4.5881 GB
*        Pipeline completed time: 2020-06-15 09:50:34
