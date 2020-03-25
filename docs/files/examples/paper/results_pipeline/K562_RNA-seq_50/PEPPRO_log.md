### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name K562_RNA-seq_50 --genome hg38 --input /project/shefflab/data/guertin/fastq/K562_50pct_RNArc_r2.fastq.gz --single-or-paired single --protocol PRO --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/dev4/results_pipeline -P 4 -M 16000`
*         Compute host:  udc-ba27-18
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/
*  Pipeline started at:   (02-27 09:23:35) elapsed: 0.0 _TIME_

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
*              `cores`:  `4`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `False`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data/guertin/fastq/K562_50pct_RNArc_r2.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `16000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed/peppro/paper/dev4/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `False`
*        `sample_name`:  `K562_RNA-seq_50`
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

Local input file: /project/shefflab/data/guertin/fastq/K562_50pct_RNArc_r2.fastq.gz

> `File_mb`	796.96	PEPPRO	_RES_

> `Read_type`	single	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-27 09:23:36) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/raw/K562_RNA-seq_50.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/K562_50pct_RNArc_r2.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/raw/K562_RNA-seq_50.fastq.gz` (184353)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 184353;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/raw/K562_RNA-seq_50.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/fastq/K562_RNA-seq_50_R1.fastq`  

> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/raw/K562_RNA-seq_50.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/fastq/K562_RNA-seq_50_R1.fastq` (184354)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 0.003GB.  
  PID: 184354;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	10000000	PEPPRO	_RES_

> `Fastq_reads`	10000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/raw/K562_RNA-seq_50.fastq.gz']

### FASTQ processing:  (02-27 09:24:07) elapsed: 32.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/fastq/K562_RNA-seq_50_R1_processed.fastq`  

> `(cutadapt -j 4 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/fastq/K562_RNA-seq_50_R1.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/fastq/K562_RNA-seq_50_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/cutadapt/K562_RNA-seq_50_R1_cutadapt.txt` (184405)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 0.16GB.  
  PID: 184405;	Command: cutadapt;	Return code: 0;	Memory used: 0.16GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/fastq/K562_RNA-seq_50_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/fastq/K562_RNA-seq_50_R1_processed.fastq` (184445,184446)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 0.16GB.  
  PID: 184445;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 184446;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/cutadapt/K562_RNA-seq_50_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	6017590.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/cutadapt/K562_RNA-seq_50_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/cutadapt/K562_RNA-seq_50_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/cutadapt/K562_RNA-seq_50_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Uninformative_adapter_reads`	124631.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	1.2463	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/fastqc/K562_RNA-seq_50_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (184478)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.16GB.  
  PID: 184478;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	9875369	PEPPRO	_RES_

> `Trim_loss_rate`	1.25	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/fastq/K562_RNA-seq_50_R1_processed.fastq` (184484)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of K562_RNA-seq_50_R1_processed.fastq
Approx 5% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 10% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 15% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 20% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 25% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 30% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 35% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 40% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 45% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 50% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 55% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 60% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 65% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 70% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 75% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 80% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 85% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 90% complete for K562_RNA-seq_50_R1_processed.fastq
Approx 95% complete for K562_RNA-seq_50_R1_processed.fastq
Analysis complete for K562_RNA-seq_50_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 0.189GB.  
  PID: 184484;	Command: fastqc;	Return code: 0;	Memory used: 0.189GB

> `FastQC report r1`	fastqc/K562_RNA-seq_50_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/fastq/processed_R1.flag` (184749)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.189GB.  
  PID: 184749;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (02-27 09:25:51) elapsed: 103.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/cutadapt/K562_RNA-seq_50_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/cutadapt/K562_RNA-seq_50_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/cutadapt` (184750)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.205GB.  
  PID: 184750;	Command: Rscript;	Return code: 0;	Memory used: 0.205GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_50_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_50_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/cutadapt/K562_RNA-seq_50_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-27 09:25:57) elapsed: 6.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/cutadapt/K562_RNA-seq_50_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/cutadapt/K562_RNA-seq_50_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/cutadapt/K562_RNA-seq_50_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/cutadapt/K562_RNA-seq_50_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/cutadapt/K562_RNA-seq_50_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2304	PEPPRO	_RES_

### Prealignments (02-27 09:25:57) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-27 09:25:57) elapsed: 0.0 _TIME_


> `(bowtie2 -p 4 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_50 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/fastq/K562_RNA-seq_50_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/prealignments/K562_RNA-seq_50_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
9875369 reads; of these:
  9875369 (100.00%) were unpaired; of these:
    9301960 (94.19%) aligned 0 times
    573409 (5.81%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
5.81% overall alignment rate

> `Aligned_reads_human_rDNA`	573409.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	5.81	PEPPRO	_RES_

### Map to genome (02-27 09:27:31) elapsed: 93.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam`  

> `bowtie2 -p 4 --very-sensitive -X 2000 --rg-id K562_RNA-seq_50 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/prealignments/K562_RNA-seq_50_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/tmpi92aywwc -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_temp.bam` (184883,184884,184885)
<pre>
9301960 reads; of these:
  9301960 (100.00%) were unpaired; of these:
    802084 (8.62%) aligned 0 times
    5741433 (61.72%) aligned exactly 1 time
    2758443 (29.65%) aligned >1 times
91.38% overall alignment rate
[bam_sort_core] merging from 2 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:13:03. Running peak memory: 3.533GB.  
  PID: 184883;	Command: bowtie2;	Return code: 0;	Memory used: 3.533GB  
  PID: 184884;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 184885;	Command: samtools;	Return code: 0;	Memory used: 0.865GB


> `samtools view -q 10 -b -@ 4 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_fail_qc.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam` (186308)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.533GB.  
  PID: 186308;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	8499876	PEPPRO	_RES_

> `QC_filtered_reads`	1395808	PEPPRO	_RES_

> `Aligned_reads`	7104068	PEPPRO	_RES_

> `Alignment_rate`	71.94	PEPPRO	_RES_

> `Total_efficiency`	71.04	PEPPRO	_RES_

> `Read_depth`	2.41	PEPPRO	_RES_

### Compress all unmapped read files (02-27 09:48:07) elapsed: 1237.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/prealignments/K562_RNA-seq_50_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/prealignments/K562_RNA-seq_50_human_rDNA_unmap.fq` (186916)
<pre>
</pre>
Command completed. Elapsed time: 0:00:50. Running peak memory: 3.533GB.  
  PID: 186916;	Command: pigz;	Return code: 0;	Memory used: 0.004GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_temp.bam` (186973)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.533GB.  
  PID: 186973;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	393154	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam` (186987)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.533GB.  
  PID: 186987;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_noMT.bam` (186994,186995,186996,186997)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.533GB.  
  PID: 186994;	Command: samtools;	Return code: 0;	Memory used: 0.003GB  
  PID: 186996;	Command: grep;	Return code: 0;	Memory used: 0.001GB  
  PID: 186995;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 186997;	Command: xargs;	Return code: 0;	Memory used: 0.026GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_noMT.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam` (187025)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 187025;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam` (187026)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.533GB.  
  PID: 187026;	Command: samtools;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-27 09:50:06) elapsed: 119.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam -c 4 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_bamQC.tsv` (187277)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/tmp_K562_RNA-seq_50_sort_dra8t2qi'
Processing with 4 cores...
Discarding 120 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrEBV']
Keeping 75 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.533GB.  
  PID: 187277;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.395GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_bamQC.tsv`

> `NRF`	0.89	PEPPRO	_RES_

> `PBC1`	0.96	PEPPRO	_RES_

> `PBC2`	36.68	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_unmap.bam`  

> `samtools view -b -@ 4 -f 4  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_unmap.bam` (187311)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.533GB.  
  PID: 187311;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools view -c -f 4 -@ 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_temp.bam`

> `Unmapped_reads`	802084	PEPPRO	_RES_

### Split BAM by strand (02-27 09:50:30) elapsed: 24.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_plus.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_plus.bam` (187331)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 3.533GB.  
  PID: 187331;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_minus.bam` (187364)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 3.533GB.  
  PID: 187364;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-27 09:51:45) elapsed: 74.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (187399)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.533GB.  
  PID: 187399;	Command: sed;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/plus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_plus_TssEnrichment.txt` (187401)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.533GB.  
  PID: 187401;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.196GB


> `TSS_coding_score`	17.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/minus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_minus_TssEnrichment.txt` (187419)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.533GB.  
  PID: 187419;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.201GB


> `TSS_non-coding_score`	4.7	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_minus_TssEnrichment.txt` (187435)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 187435;	Command: Rscript;	Return code: 0;	Memory used: 0.234GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_50_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_50_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt` (187453,187454,187455,187456)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 187453;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 187455;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 187454;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 187456;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_keep.txt` (187458)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 187458;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-27 09:52:08) elapsed: 23.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/hg38_ensembl_tss.bed` (187460,187461)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.533GB.  
  PID: 187460;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 187461;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/hg38_ensembl_gene_body.bed` (187466,187467)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.533GB.  
  PID: 187467;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 187466;	Command: grep;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_TSS_density.bed` (187470,187471,187472,187473)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.533GB.  
  PID: 187473;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 187470;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 187472;	Command: sort;	Return code: 0;	Memory used: 0.007GB  
  PID: 187471;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_gene_body_density.bed` (187490,187491,187492)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.533GB.  
  PID: 187490;	Command: bedtools;	Return code: 0;	Memory used: 0.042GB  
  PID: 187492;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 187491;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_TSS_density.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_pause_index.bed` (187509,187510,187511)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 187509;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 187511;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 187510;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	12.83	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_pause_index.bed` (187516)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.533GB.  
  PID: 187516;	Command: Rscript;	Return code: 0;	Memory used: 0.235GB

> `Pause index`	QC_hg38/K562_RNA-seq_50_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_50_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_pause_index.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_pause_index.bed` (187533)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 187533;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (02-27 09:52:51) elapsed: 43.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_plus.bam`
7104068 2678145

> `Plus_FRiP`	0.38	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_minus.bam`
7104068 2628776

> `Minus_FRiP`	0.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/signal_hg38/K562_RNA-seq_50_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/hg38_gene_sort.bed` (187562,187563)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 187563;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 187562;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/signal_hg38/K562_RNA-seq_50_gene_coverage.bed` (187566)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.533GB.  
  PID: 187566;	Command: bedtools;	Return code: 0;	Memory used: 0.042GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/raw/hg38_annotations.bed.gz` (187580)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 187580;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/raw/hg38_annotations.bed` (187581)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 187581;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (02-27 09:53:21) elapsed: 31.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/raw/hg38_annotations.bed` (187595)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.533GB.  
  PID: 187595;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/3_UTR"` (187598)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 187598;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/3_UTR_sort.bed` (187599,187600,187601,187602)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 187599;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 187600;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 187602;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 187601;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_3_UTR_plus_coverage.bed` (187604)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 187604;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_3_UTR_minus_coverage.bed` (187612)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 187612;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/5_UTR"` (187619)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 187619;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/5_UTR_sort.bed` (187620,187621,187622,187623)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 187620;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 187621;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 187623;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 187622;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_5_UTR_plus_coverage.bed` (187625)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 187625;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_5_UTR_minus_coverage.bed` (187632)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 187632;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Enhancer_sort.bed` (187639,187640,187641,187642)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 187639;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 187640;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 187642;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB  
  PID: 187641;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Enhancer_plus_coverage.bed` (187645)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 187645;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Enhancer_minus_coverage.bed` (187653)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 187653;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Exon_sort.bed` (187660,187661,187662,187663)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.533GB.  
  PID: 187660;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 187661;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 187663;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 187662;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Exon_plus_coverage.bed` (187668)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.533GB.  
  PID: 187668;	Command: bedtools;	Return code: 0;	Memory used: 0.025GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Exon_minus_coverage.bed` (187676)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.533GB.  
  PID: 187676;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Intron_sort.bed` (187687,187688,187689,187690)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.533GB.  
  PID: 187687;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 187689;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 187688;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 187690;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Intron_plus_coverage.bed` (187693)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.533GB.  
  PID: 187693;	Command: bedtools;	Return code: 0;	Memory used: 0.027GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Intron_minus_coverage.bed` (187702)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.533GB.  
  PID: 187702;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter_sort.bed` (187710,187711,187712,187713)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 187710;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 187711;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 187713;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 187712;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Promoter_plus_coverage.bed` (187716)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 187716;	Command: bedtools;	Return code: 0;	Memory used: 0.023GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Promoter_minus_coverage.bed` (187723)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 187723;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter_Flanking_Region"` (187730)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 187730;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter_Flanking_Region_sort.bed` (187731,187732,187733,187734)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 187731;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 187733;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 187732;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 187734;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Promoter_Flanking_Region_plus_coverage.bed` (187738)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.533GB.  
  PID: 187738;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Promoter_Flanking_Region_minus_coverage.bed` (187936)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.533GB.  
  PID: 187936;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB


### Plot cFRiF/FRiF (02-27 09:55:11) elapsed: 109.0 _TIME_


> `samtools view -@ 4 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_50 -z 3099922541 -n 3418574 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Promoter_Flanking_Region_plus_coverage.bed` (187955)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.533GB.  
  PID: 187955;	Command: Rscript;	Return code: 0;	Memory used: 0.529GB

> `cFRiF`	QC_hg38/K562_RNA-seq_50_cFRiF.pdf	cFRiF	QC_hg38/K562_RNA-seq_50_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_50 -z 3099922541 -n 3418574 -y frif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_Promoter_Flanking_Region_plus_coverage.bed` (187999)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.533GB.  
  PID: 187999;	Command: Rscript;	Return code: 0;	Memory used: 0.529GB

> `FRiF`	QC_hg38/K562_RNA-seq_50_FRiF.pdf	FRiF	QC_hg38/K562_RNA-seq_50_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-27 09:56:14) elapsed: 63.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/hg38_exons_sort.bed` (188032,188033)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.533GB.  
  PID: 188033;	Command: bedtools;	Return code: 0;	Memory used: 0.08GB  
  PID: 188032;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/hg38_introns_sort.bed` (188041,188042,188043)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.533GB.  
  PID: 188042;	Command: bedtools;	Return code: 0;	Memory used: 0.078GB  
  PID: 188041;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 188043;	Command: bedtools;	Return code: 0;	Memory used: 0.092GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_exons_coverage.bed` (188055)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.533GB.  
  PID: 188055;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_introns_coverage.bed` (188068)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.533GB.  
  PID: 188068;	Command: bedtools;	Return code: 0;	Memory used: 0.032GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.104068)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_exons_rpkm.bed` (188081,188082,188083)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 188081;	Command: awk;	Return code: 0;	Memory used: 0.005GB  
  PID: 188083;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 188082;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.104068)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_introns_rpkm.bed` (188086,188087,188088)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.533GB.  
  PID: 188086;	Command: awk;	Return code: 0;	Memory used: 0.006GB  
  PID: 188088;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 188087;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_introns_rpkm.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_exon_intron_ratios.bed` (188090,188091,188092)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 188090;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 188092;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 188091;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	4.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_exon_intron_ratios.bed --annotate` (188098)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.533GB.  
  PID: 188098;	Command: Rscript;	Return code: 0;	Memory used: 0.235GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_50_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_50_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_exon_intron_ratios.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/QC_hg38/K562_RNA-seq_50_exon_intron_ratios.bed` (188116)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.533GB.  
  PID: 188116;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Produce bigWig files (02-27 09:57:03) elapsed: 49.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/signal_hg38/K562_RNA-seq_50_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/signal_hg38/K562_RNA-seq_50_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_plus.bam` (188123)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.533GB.  
  PID: 188123;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/signal_hg38/K562_RNA-seq_50_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/signal_hg38/K562_RNA-seq_50_plus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge` (188127)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_50_plus_cuttrace_px1xi48z'
Processing with 1 cores...
Discarding 130 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_GL000218v1', 'chrEBV']
Keeping 65 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270582v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2']
Reduce step (merge files)...
Merging 65 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/signal_hg38/K562_RNA-seq_50_plus_exact_body_0-mer.bw'
Merging 65 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/signal_hg38/K562_RNA-seq_50_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:16:08. Running peak memory: 3.533GB.  
  PID: 188127;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.484GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/signal_hg38/K562_RNA-seq_50_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/signal_hg38/K562_RNA-seq_50_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_minus.bam` (190259)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.533GB.  
  PID: 190259;	Command: samtools;	Return code: 0;	Memory used: 0.005GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/signal_hg38/K562_RNA-seq_50_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/signal_hg38/K562_RNA-seq_50_minus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge` (190263)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/aligned_hg38/K562_RNA-seq_50_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_50_minus_cuttrace_1tne8kw3'
Processing with 1 cores...
stdin is empty of data
Discarding 128 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr16_KI270728v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 67 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr15_KI270727v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 67 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/signal_hg38/K562_RNA-seq_50_minus_exact_body_0-mer.bw'
Merging 67 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_50/signal_hg38/K562_RNA-seq_50_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:16:01. Running peak memory: 3.533GB.  
  PID: 190263;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.481GB

Starting cleanup: 52 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  1:05:47
*  Total elapsed time (all runs):  1:22:09
*         Peak memory (this run):  3.5329 GB
*        Pipeline completed time: 2020-02-27 10:29:22
