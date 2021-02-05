### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_50 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_50pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 16 -M 16000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aw29-26b
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/
*  Pipeline started at:   (06-15 07:11:05) elapsed: 3.0 _TIME_

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
*              `cores`:  `16`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_50pct.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `16000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `K562_PRO-seq_50`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_50pct.fastq.gz

> `File_mb`	18213.9	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:11:05) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/raw/K562_PRO-seq_50.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_PRO_50pct.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/raw/K562_PRO-seq_50.fastq.gz` (208284)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 208284;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/raw/K562_PRO-seq_50.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastq/K562_PRO-seq_50_R1.fastq`  

> `pigz -f -p 16 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/raw/K562_PRO-seq_50.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastq/K562_PRO-seq_50_R1.fastq` (208287)
<pre>
</pre>
Command completed. Elapsed time: 0:09:50. Running peak memory: 0.003GB.  
  PID: 208287;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	248316549	PEPPRO	_RES_

> `Fastq_reads`	248316549	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/raw/K562_PRO-seq_50.fastq.gz']

### FASTQ processing:  (06-15 07:32:58) elapsed: 1313.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastq/K562_PRO-seq_50_R1_processed.fastq`  

> `(cutadapt -j 16 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastq/K562_PRO-seq_50_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastq/K562_PRO-seq_50_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastq/K562_PRO-seq_50_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/cutadapt/K562_PRO-seq_50_R1_cutadapt.txt` (210371)
<pre>
</pre>
Command completed. Elapsed time: 0:11:40. Running peak memory: 6.72GB.  
  PID: 210371;	Command: cutadapt;	Return code: 0;	Memory used: 6.72GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastq/K562_PRO-seq_50_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastq/K562_PRO-seq_50_R1_processed.fastq` (211485,211486)
<pre>
</pre>
Command completed. Elapsed time: 0:06:14. Running peak memory: 6.72GB.  
  PID: 211486;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB  
  PID: 211485;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/cutadapt/K562_PRO-seq_50_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	211550349.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/cutadapt/K562_PRO-seq_50_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/cutadapt/K562_PRO-seq_50_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastq/K562_PRO-seq_50_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	6217636.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	2.5039	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastqc/K562_PRO-seq_50_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (212485)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.72GB.  
  PID: 212485;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	242098913	PEPPRO	_RES_

> `Trim_loss_rate`	2.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastq/K562_PRO-seq_50_R1_processed.fastq` (212735)
<pre>
Started analysis of K562_PRO-seq_50_R1_processed.fastq
Approx 5% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 10% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 15% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 20% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 25% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 30% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 35% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 40% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 45% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 50% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 55% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 60% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 65% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 70% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 75% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 80% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 85% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 90% complete for K562_PRO-seq_50_R1_processed.fastq
Approx 95% complete for K562_PRO-seq_50_R1_processed.fastq
Analysis complete for K562_PRO-seq_50_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:13:13. Running peak memory: 6.72GB.  
  PID: 212735;	Command: fastqc;	Return code: 0;	Memory used: 0.252GB

> `FastQC report r1`	fastqc/K562_PRO-seq_50_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastq/processed_R1.flag` (213985)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.72GB.  
  PID: 213985;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-15 08:09:08) elapsed: 2170.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/cutadapt/K562_PRO-seq_50_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/cutadapt/K562_PRO-seq_50_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/cutadapt` (213991)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.72GB.  
  PID: 213991;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_50_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_50_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/cutadapt/K562_PRO-seq_50_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 08:09:13) elapsed: 5.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/cutadapt/K562_PRO-seq_50_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/cutadapt/K562_PRO-seq_50_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/cutadapt/K562_PRO-seq_50_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/cutadapt/K562_PRO-seq_50_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/cutadapt/K562_PRO-seq_50_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2313	PEPPRO	_RES_

### Prealignments (06-15 08:09:13) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 08:09:13) elapsed: 0.0 _TIME_


> `(bowtie2 -p 16 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_PRO-seq_50 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/fastq/K562_PRO-seq_50_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/prealignments/K562_PRO-seq_50_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
242098913 reads; of these:
  242098913 (100.00%) were unpaired; of these:
    219756437 (90.77%) aligned 0 times
    22342476 (9.23%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.23% overall alignment rate

> `Aligned_reads_human_rDNA`	22342476.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.23	PEPPRO	_RES_

### Map to genome (06-15 08:45:23) elapsed: 2169.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam`  

> `bowtie2 -p 16 --very-sensitive --rg-id K562_PRO-seq_50 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/prealignments/K562_PRO-seq_50_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/tmph50h2hb_ -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_temp.bam` (219931,220017,220018)
<pre>
219756437 reads; of these:
  219756437 (100.00%) were unpaired; of these:
    2697646 (1.23%) aligned 0 times
    164092013 (74.67%) aligned exactly 1 time
    52966778 (24.10%) aligned >1 times
98.77% overall alignment rate
[bam_sort_core] merging from 70 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 2:07:48. Running peak memory: 6.72GB.  
  PID: 220017;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 219931;	Command: bowtie2;	Return code: 0;	Memory used: 3.845GB  
  PID: 220018;	Command: samtools;	Return code: 0;	Memory used: 0.9GB


> `samtools view -q 10 -b -@ 16 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam` (237040)
<pre>
</pre>
Command completed. Elapsed time: 0:06:04. Running peak memory: 6.72GB.  
  PID: 237040;	Command: samtools;	Return code: 0;	Memory used: 0.026GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	217058791	PEPPRO	_RES_

> `QC_filtered_reads`	23375045	PEPPRO	_RES_

> `Aligned_reads`	193683746	PEPPRO	_RES_

> `Alignment_rate`	80.0	PEPPRO	_RES_

> `Total_efficiency`	78.0	PEPPRO	_RES_

> `Read_depth`	16.03	PEPPRO	_RES_

### Compress all unmapped read files (06-15 11:39:29) elapsed: 10447.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/prealignments/K562_PRO-seq_50_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/prealignments/K562_PRO-seq_50_human_rDNA_unmap.fq` (241954)
<pre>
</pre>
Command completed. Elapsed time: 0:05:30. Running peak memory: 6.72GB.  
  PID: 241954;	Command: pigz;	Return code: 0;	Memory used: 0.014GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_temp.bam` (242600)
<pre>
</pre>
Command completed. Elapsed time: 0:03:22. Running peak memory: 6.72GB.  
  PID: 242600;	Command: samtools;	Return code: 0;	Memory used: 0.022GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	4576683	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam` (243016)
<pre>
</pre>
Command completed. Elapsed time: 0:02:56. Running peak memory: 6.72GB.  
  PID: 243016;	Command: samtools;	Return code: 0;	Memory used: 0.022GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/chr_sizes.bed` (270996,270997,270998,270999)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.72GB.  
  PID: 270997;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 270999;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 270996;	Command: samtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 270998;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/chr_sizes.bed -b -@ 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_noMT.bam` (271001)
<pre>
</pre>
Command completed. Elapsed time: 0:03:37. Running peak memory: 6.72GB.  
  PID: 271001;	Command: samtools;	Return code: 0;	Memory used: 0.023GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam` (279626)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.72GB.  
  PID: 279626;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam` (280008)
<pre>
</pre>
Command completed. Elapsed time: 0:02:56. Running peak memory: 6.72GB.  
  PID: 280008;	Command: samtools;	Return code: 0;	Memory used: 0.022GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 12:06:45) elapsed: 1635.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam -c 16 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_bamQC.tsv` (338745)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/tmp_K562_PRO-seq_50_sort__n9_22y7'
Processing with 16 cores...
Discarding 86 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 109 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:03:15. Running peak memory: 7.423GB.  
  PID: 338745;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 7.423GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_bamQC.tsv`

> `NRF`	0.47	PEPPRO	_RES_

> `PBC1`	0.7	PEPPRO	_RES_

> `PBC2`	4.99	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_unmap.bam`  

> `samtools view -b -@ 16 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_unmap.bam` (382983)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 7.423GB.  
  PID: 382983;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools view -c -f 4 -@ 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_temp.bam`

> `Unmapped_reads`	2697646	PEPPRO	_RES_

### Split BAM by strand (06-15 12:11:11) elapsed: 267.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_plus.bam` (383351)
<pre>
</pre>
Command completed. Elapsed time: 0:14:01. Running peak memory: 7.423GB.  
  PID: 383351;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_minus.bam` (51589)
<pre>
</pre>
Command completed. Elapsed time: 0:13:29. Running peak memory: 7.423GB.  
  PID: 51589;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 12:38:41) elapsed: 1649.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (172739)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 172739;	Command: sed;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/plus_TSS.tsv -p ends -c 16 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_plus_TssEnrichment.txt` (173003)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 7.423GB.  
  PID: 173003;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.063GB


> `TSS_coding_score`	13.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/minus_TSS.tsv -p ends -c 16 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_minus_TssEnrichment.txt` (181856)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 7.423GB.  
  PID: 181856;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.049GB


> `TSS_non-coding_score`	5.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_minus_TssEnrichment.txt` (187364)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 7.423GB.  
  PID: 187364;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `TSS enrichment`	QC_hg38/K562_PRO-seq_50_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_50_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt` (187419,187420,187421,187422)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 187419;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 187421;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 187420;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 187422;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_keep.txt` (187424)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 187424;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 12:39:33) elapsed: 53.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/hg38_ensembl_tss.bed` (187426,187427)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 7.423GB.  
  PID: 187426;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 187427;	Command: bedtools;	Return code: 0;	Memory used: 0.095GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/hg38_ensembl_gene_body.bed` (187432,187433)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 187432;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 187433;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_TSS_density.bed` (187435,187436,187437,187438)
<pre>
</pre>
Command completed. Elapsed time: 0:04:10. Running peak memory: 7.423GB.  
  PID: 187435;	Command: bedtools;	Return code: 0;	Memory used: 0.055GB  
  PID: 187437;	Command: sort;	Return code: 0;	Memory used: 0.014GB  
  PID: 187436;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 187438;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_gene_body_density.bed` (194172,194173,194175)
<pre>
</pre>
Command completed. Elapsed time: 0:06:38. Running peak memory: 7.423GB.  
  PID: 194175;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 194172;	Command: bedtools;	Return code: 0;	Memory used: 0.524GB  
  PID: 194173;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/tmpv1edhois` (195123,195124,195125)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 195123;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 195125;	Command: env;	Return code: 0;	Memory used: 0.005GB  
  PID: 195124;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/tmpv1edhois | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.119401) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/tmpv1edhois > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_pause_index.bed` (195131)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 195131;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.88	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_pause_index.bed` (195136)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 7.423GB.  
  PID: 195136;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Pause index`	QC_hg38/K562_PRO-seq_50_pause_index.pdf	Pause index	QC_hg38/K562_PRO-seq_50_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_pause_index.bed.gz`  

> `pigz -f -p 16 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_pause_index.bed` (195160)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 195160;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 12:50:29) elapsed: 656.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_plus.bam`
193683746 69078220

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_minus.bam`
193683746 65709582

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/signal_hg38/K562_PRO-seq_50_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/hg38_gene_sort.bed` (195667,195668)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.423GB.  
  PID: 195667;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 195668;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/signal_hg38/K562_PRO-seq_50_gene_coverage.bed` (195670)
<pre>
</pre>
Command completed. Elapsed time: 0:06:33. Running peak memory: 7.423GB.  
  PID: 195670;	Command: bedtools;	Return code: 0;	Memory used: 0.52GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/raw/hg38_annotations.bed.gz` (196281)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 196281;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 16 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/raw/hg38_annotations.bed` (196287)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 196287;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 13:02:37) elapsed: 728.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/raw/hg38_annotations.bed` (196296)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 7.423GB.  
  PID: 196296;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Enhancer_sort.bed` (196300,196301,196302,196303)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.423GB.  
  PID: 196300;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 196301;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 196303;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 196302;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Enhancer_plus_coverage.bed` (196306)
<pre>
</pre>
Command completed. Elapsed time: 0:02:07. Running peak memory: 7.423GB.  
  PID: 196306;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Enhancer_minus_coverage.bed` (196423)
<pre>
</pre>
Command completed. Elapsed time: 0:02:00. Running peak memory: 7.423GB.  
  PID: 196423;	Command: bedtools;	Return code: 0;	Memory used: 0.038GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter_sort.bed` (196718,196719,196720,196721)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 196718;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 196719;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 196721;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 196720;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Promoter_plus_coverage.bed` (196723)
<pre>
</pre>
Command completed. Elapsed time: 0:02:12. Running peak memory: 7.423GB.  
  PID: 196723;	Command: bedtools;	Return code: 0;	Memory used: 0.284GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Promoter_minus_coverage.bed` (196875)
<pre>
</pre>
Command completed. Elapsed time: 0:02:10. Running peak memory: 7.423GB.  
  PID: 196875;	Command: bedtools;	Return code: 0;	Memory used: 0.239GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter_Flanking_Region"` (197225)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 197225;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter_Flanking_Region_sort.bed` (197226,197227,197228,197229)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.423GB.  
  PID: 197226;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 197228;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 197227;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 197229;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Promoter_Flanking_Region_plus_coverage.bed` (197232)
<pre>
</pre>
Command completed. Elapsed time: 0:02:06. Running peak memory: 7.423GB.  
  PID: 197232;	Command: bedtools;	Return code: 0;	Memory used: 0.053GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Promoter_Flanking_Region_minus_coverage.bed` (199936)
<pre>
</pre>
Command completed. Elapsed time: 0:02:07. Running peak memory: 7.423GB.  
  PID: 199936;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/5_UTR"` (201131)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 201131;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/5_UTR_sort.bed` (201132,201133,201134,201135)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.423GB.  
  PID: 201132;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 201133;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 201135;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 201134;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_5_UTR_plus_coverage.bed` (201138)
<pre>
</pre>
Command completed. Elapsed time: 0:02:03. Running peak memory: 7.423GB.  
  PID: 201138;	Command: bedtools;	Return code: 0;	Memory used: 0.042GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_5_UTR_minus_coverage.bed` (201246)
<pre>
</pre>
Command completed. Elapsed time: 0:01:59. Running peak memory: 7.423GB.  
  PID: 201246;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/3_UTR"` (201351)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 201351;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/3_UTR_sort.bed` (201352,201353,201354,201355)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.423GB.  
  PID: 201352;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 201353;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 201355;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 201354;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_3_UTR_plus_coverage.bed` (201358)
<pre>
</pre>
Command completed. Elapsed time: 0:02:07. Running peak memory: 7.423GB.  
  PID: 201358;	Command: bedtools;	Return code: 0;	Memory used: 0.056GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_3_UTR_minus_coverage.bed` (201699)
<pre>
</pre>
Command completed. Elapsed time: 0:02:02. Running peak memory: 7.423GB.  
  PID: 201699;	Command: bedtools;	Return code: 0;	Memory used: 0.084GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Exon_sort.bed` (201806,201807,201808,201809)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 7.423GB.  
  PID: 201806;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 201807;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 201809;	Command: bedtools;	Return code: 0;	Memory used: 0.174GB  
  PID: 201808;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Exon_plus_coverage.bed` (201814)
<pre>
</pre>
Command completed. Elapsed time: 0:02:16. Running peak memory: 7.423GB.  
  PID: 201814;	Command: bedtools;	Return code: 0;	Memory used: 0.239GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Exon_minus_coverage.bed` (202129)
<pre>
</pre>
Command completed. Elapsed time: 0:02:07. Running peak memory: 7.423GB.  
  PID: 202129;	Command: bedtools;	Return code: 0;	Memory used: 0.097GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Intron_sort.bed` (202347,202348,202349,202350)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.423GB.  
  PID: 202347;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 202349;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 202348;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 202350;	Command: bedtools;	Return code: 0;	Memory used: 0.083GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Intron_plus_coverage.bed` (202361)
<pre>
</pre>
Command completed. Elapsed time: 0:02:33. Running peak memory: 7.423GB.  
  PID: 202361;	Command: bedtools;	Return code: 0;	Memory used: 0.228GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Intron_minus_coverage.bed` (202937)
<pre>
</pre>
Command completed. Elapsed time: 0:02:29. Running peak memory: 7.423GB.  
  PID: 202937;	Command: bedtools;	Return code: 0;	Memory used: 0.222GB


### Plot cFRiF/FRiF (06-15 13:33:06) elapsed: 1830.0 _TIME_


> `samtools view -@ 16 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_50 -z 3099922541 -n 96376229 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Intron_plus_coverage.bed` (203106)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 7.423GB.  
  PID: 203106;	Command: Rscript;	Return code: 0;	Memory used: 0.461GB

> `cFRiF`	QC_hg38/K562_PRO-seq_50_cFRiF.pdf	cFRiF	QC_hg38/K562_PRO-seq_50_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_50 -z 3099922541 -n 96376229 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_Intron_plus_coverage.bed` (203154)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 7.423GB.  
  PID: 203154;	Command: Rscript;	Return code: 0;	Memory used: 0.467GB

> `FRiF`	QC_hg38/K562_PRO-seq_50_FRiF.pdf	FRiF	QC_hg38/K562_PRO-seq_50_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 13:34:26) elapsed: 80.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/hg38_exons_sort.bed` (203188,203189)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 7.423GB.  
  PID: 203189;	Command: bedtools;	Return code: 0;	Memory used: 0.098GB  
  PID: 203188;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/hg38_introns_sort.bed` (203195,203196,203197)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 7.423GB.  
  PID: 203195;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 203197;	Command: bedtools;	Return code: 0;	Memory used: 0.003GB  
  PID: 203196;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_exons_coverage.bed` (203203)
<pre>
</pre>
Command completed. Elapsed time: 0:04:27. Running peak memory: 7.423GB.  
  PID: 203203;	Command: bedtools;	Return code: 0;	Memory used: 0.15GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_introns_coverage.bed` (203723)
<pre>
</pre>
Command completed. Elapsed time: 0:05:20. Running peak memory: 7.423GB.  
  PID: 203723;	Command: bedtools;	Return code: 0;	Memory used: 0.251GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/193.683746)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_exons_rpkm.bed` (204334,204335,204336)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.423GB.  
  PID: 204334;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 204336;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 204335;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/193.683746)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_introns_rpkm.bed` (204338,204339,204340)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.423GB.  
  PID: 204338;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 204340;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 204339;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_exon_intron_ratios.bed` (204343,204344,204345)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 204343;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 204345;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 204344;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_exon_intron_ratios.bed --annotate` (204351)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 7.423GB.  
  PID: 204351;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/K562_PRO-seq_50_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_PRO-seq_50_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_exon_intron_ratios.bed.gz`  

> `pigz -f -p 16 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/QC_hg38/K562_PRO-seq_50_exon_intron_ratios.bed` (204374)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.423GB.  
  PID: 204374;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 13:44:30) elapsed: 604.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/signal_hg38/K562_PRO-seq_50_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/signal_hg38/K562_PRO-seq_50_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_plus.bam` (204383)
<pre>
</pre>
Command completed. Elapsed time: 0:01:22. Running peak memory: 7.423GB.  
  PID: 204383;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/signal_hg38/K562_PRO-seq_50_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/signal_hg38/K562_PRO-seq_50_plus_smooth_body_0-mer.bw -p 10 --variable-step --tail-edge --scale 193683746.0` (204651)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_plus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j'
Processing with 5 cores...
stdin is empty of data
stdin is empty of data
mustOpen: Can't open tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr11.txt.bw to write: No such file or directory
mustOpen: Can't open tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr3.txt_smooth.bw to write: No such file or directory
tee: tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr4.txt_cuts.txt: No such file or directory
mustOpen: Can't open tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr11.txt_smooth.bw to write: No such file or directory
tee: tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr12.txt_cuts.txt: No such file or directory
mustOpen: Can't open tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr4.txt.bw to write: No such file or directory
mustOpen: Can't open tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr12.txt.bw to write: No such file or directory
mustOpen: Can't open tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr12.txt_smooth.bw to write: No such file or directory
mustOpen: Can't open tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr4.txt_smooth.bw to write: No such file or directory
tee: tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr5.txt_cuts.txt: No such file or directory
mustOpen: Can't open tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr5.txt.bw to write: No such file or directory
mustOpen: Can't open tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr5.txt_smooth.bw to write: No such file or directory
tee: tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr6.txt_cuts.txt: No such file or directory
mustOpen: Can't open tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr6.txt.bw to write: No such file or directory
mustOpen: Can't open tmp_K562_PRO-seq_50_plus_cuttrace_hliil13j/chr6.txt_smooth.bw to write: No such file or directory
Discarding 92 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1']
Keeping 103 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 103 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/signal_hg38/K562_PRO-seq_50_plus_exact_body_0-mer.bw'
bigWigCat v 4 - merge non-overlapping bigWig files
directly into bigWig format
usage:
   bigWigCat out.bw in1.bw in2.bw ...
Where in*.bw is in big wig format
and out.bw is the output indexed big wig file.
options:
   -itemsPerSlot=N - Number of data points bundled at lowest level. Default 1024

Note: must use wigToBigWig -fixedSummaries -keepAllChromosomes (perhaps in parallel cluster jobs) to create the input files.
Note: By non-overlapping we mean the entire span of each file, from first data point to last data point, must not overlap with that of other files.

Merging 103 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/signal_hg38/K562_PRO-seq_50_plus_smooth_body_0-mer.bw'
bigWigCat v 4 - merge non-overlapping bigWig files
directly into bigWig format
usage:
   bigWigCat out.bw in1.bw in2.bw ...
Where in*.bw is in big wig format
and out.bw is the output indexed big wig file.
options:
   -itemsPerSlot=N - Number of data points bundled at lowest level. Default 1024

Note: must use wigToBigWig -fixedSummaries -keepAllChromosomes (perhaps in parallel cluster jobs) to create the input files.
Note: By non-overlapping we mean the entire span of each file, from first data point to last data point, must not overlap with that of other files.

</pre>
Command completed. Elapsed time: 0:07:57. Running peak memory: 7.423GB.  
  PID: 204651;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.644GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/signal_hg38/K562_PRO-seq_50_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/signal_hg38/K562_PRO-seq_50_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_minus.bam` (206359)
<pre>
</pre>
Command completed. Elapsed time: 0:01:25. Running peak memory: 7.423GB.  
  PID: 206359;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/signal_hg38/K562_PRO-seq_50_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/signal_hg38/K562_PRO-seq_50_minus_smooth_body_0-mer.bw -p 10 --variable-step --tail-edge --scale 193683746.0` (206630)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/aligned_hg38/K562_PRO-seq_50_minus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_50_minus_cuttrace_6qrk7or_'
Processing with 5 cores...
stdin is empty of data
stdin is empty of data
Discarding 98 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 97 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 97 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/signal_hg38/K562_PRO-seq_50_minus_exact_body_0-mer.bw'
Merging 97 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_50/signal_hg38/K562_PRO-seq_50_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:49. Running peak memory: 7.423GB.  
  PID: 206630;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.583GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  6:53:02
*  Total elapsed time (all runs):  8:51:51
*         Peak memory (this run):  7.4232 GB
*        Pipeline completed time: 2020-06-15 14:04:04
