### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name Jurkat_ChRO-seq_1 --genome hg38 --input /project/shefflab/data/sra_fastq/SRR7616133.fastq.gz --single-or-paired single --protocol PRO --umi-len 6 --max-len -1 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-aj38-13c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/
*  Pipeline started at:   (02-18 11:08:58) elapsed: 0.0 _TIME_

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
*              `cores`:  `8`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `False`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data/sra_fastq/SRR7616133.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `16000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed/peppro/paper/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `False`
*        `sample_name`:  `Jurkat_ChRO-seq_1`
*              `scale`:  `False`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `single`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `6`
*          `verbosity`:  `None`

----------------------------------------

Local input file: /project/shefflab/data/sra_fastq/SRR7616133.fastq.gz

> `File_mb`	1330.25	PEPPRO	_RES_

> `Read_type`	single	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-18 11:08:58) elapsed: 0.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz`  

> `ln -sf /project/shefflab/data/sra_fastq/SRR7616133.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz` (327686)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 327686;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1.fastq` (327689)
<pre>
</pre>
Command completed. Elapsed time: 0:02:30. Running peak memory: 0.003GB.  
  PID: 327689;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	37740511	PEPPRO	_RES_

> `Fastq_reads`	37740511	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz']

### FASTQ processing:  (02-18 11:11:52) elapsed: 173.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_processed.fastq`  

> `(cutadapt -j 8 -m 8 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt` (328109)
<pre>
</pre>
Command completed. Elapsed time: 0:01:18. Running peak memory: 0.289GB.  
  PID: 328109;	Command: cutadapt;	Return code: 0;	Memory used: 0.289GB


> `seqtk trimfq -b 6 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_noadap.fastq | seqtk seq -L 8 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_processed.fastq` (328231,328232)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 0.289GB.  
  PID: 328231;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 328232;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	32783879	PEPPRO	_RES_

> `Trim_loss_rate`	13.13	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_processed.fastq` (328289)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 5% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 10% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 15% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 20% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 25% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 30% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 35% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 40% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 45% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 50% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 55% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 60% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 65% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 70% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 75% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 80% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 85% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 90% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 95% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Analysis complete for Jurkat_ChRO-seq_1_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:49. Running peak memory: 0.289GB.  
  PID: 328289;	Command: fastqc;	Return code: 0;	Memory used: 0.184GB

> `FastQC report r1`	fastqc/Jurkat_ChRO-seq_1_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_dedup.fastq /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_noadap.fastq` (329801)
<pre>
[INFO][0m 14846147 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 2.009GB.  
  PID: 329801;	Command: seqkit;	Return code: 0;	Memory used: 2.009GB


> `seqtk trimfq -b 6 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_dedup.fastq | seqtk seq -L 8 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_trimmed.fastq` (329859,329860)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 2.009GB.  
  PID: 329859;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 329860;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	34473620.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_too_short`	2012355.0	PEPPRO	_RES_

> `Duplicate_reads`	14846147.0	PEPPRO	_RES_

> `Pct_reads_too_short`	5.3321	PEPPRO	_RES_

### Plot adapter insertion distribution (02-18 11:17:20) elapsed: 328.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/cutadapt -u 6` (329927)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 2.009GB.  
  PID: 329927;	Command: Rscript;	Return code: 0;	Memory used: 0.203GB

> `Adapter insertion distribution`	cutadapt/Jurkat_ChRO-seq_1_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/Jurkat_ChRO-seq_1_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-18 11:17:26) elapsed: 6.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.9566	PEPPRO	_RES_

### Prealignments (02-18 11:17:26) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-18 11:17:26) elapsed: 0.0 _TIME_


> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id Jurkat_ChRO-seq_1 -U /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
32783879 reads; of these:
  32783879 (100.00%) were unpaired; of these:
    29290636 (89.34%) aligned 0 times
    3493243 (10.66%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.66% overall alignment rate

> `Aligned_reads_human_rDNA`	3493243.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.66	PEPPRO	_RES_

### Map to human_rDNA (02-18 11:21:32) elapsed: 245.0 _TIME_


> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id Jurkat_ChRO-seq_1 -U /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap_dups.fq 2>&1 > /dev/null)`

### Map to genome (02-18 11:23:31) elapsed: 119.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id Jurkat_ChRO-seq_1 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/tmpellgvmi2 -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp.bam` (330654,330655,330656)
<pre>
29290636 reads; of these:
  29290636 (100.00%) were unpaired; of these:
    1113286 (3.80%) aligned 0 times
    17309978 (59.10%) aligned exactly 1 time
    10867372 (37.10%) aligned >1 times
96.20% overall alignment rate
[bam_sort_core] merging from 7 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:17:52. Running peak memory: 3.606GB.  
  PID: 330655;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 330654;	Command: bowtie2;	Return code: 0;	Memory used: 3.606GB  
  PID: 330656;	Command: samtools;	Return code: 0;	Memory used: 0.881GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam` (332661)
<pre>
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 3.606GB.  
  PID: 332661;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	28177350	PEPPRO	_RES_

> `QC_filtered_reads`	6839271	PEPPRO	_RES_

> `Aligned_reads`	21338079	PEPPRO	_RES_

> `Alignment_rate`	65.09	PEPPRO	_RES_

> `Total_efficiency`	56.54	PEPPRO	_RES_

> `Read_depth`	3.87	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id Jurkat_ChRO-seq_1 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap_dups.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/tmpellgvmi2 -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp_dups.bam` (333755,333760,333761)
<pre>
18024676 reads; of these:
  18024676 (100.00%) were unpaired; of these:
    856458 (4.75%) aligned 0 times
    10720803 (59.48%) aligned exactly 1 time
    6447415 (35.77%) aligned >1 times
95.25% overall alignment rate
[bam_sort_core] merging from 4 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:10:52. Running peak memory: 3.606GB.  
  PID: 333760;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 333755;	Command: bowtie2;	Return code: 0;	Memory used: 3.605GB  
  PID: 333761;	Command: samtools;	Return code: 0;	Memory used: 0.875GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp_dups.bam > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort_dups.bam` (334850)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 3.606GB.  
  PID: 334850;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


### Compress all unmapped read files (02-18 12:02:03) elapsed: 2312.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap.fq` (334891)
<pre>
</pre>
Command completed. Elapsed time: 0:00:43. Running peak memory: 3.606GB.  
  PID: 334891;	Command: pigz;	Return code: 0;	Memory used: 0.008GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp.bam` (334946)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.606GB.  
  PID: 334946;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	4738	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam` (334970)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.606GB.  
  PID: 334970;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_noMT.bam` (334983,334984,334985,334986)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.606GB.  
  PID: 334984;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 334983;	Command: samtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 334985;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 334986;	Command: xargs;	Return code: 0;	Memory used: 0.047GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam` (335013)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.606GB.  
  PID: 335013;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam` (335015)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.606GB.  
  PID: 335015;	Command: samtools;	Return code: 0;	Memory used: 0.015GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	70	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp_dups.bam` (335071)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.606GB.  
  PID: 335071;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort_dups.bam`  

### Calculate library complexity (02-18 12:04:47) elapsed: 164.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_out.txt -B /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort_dups.bam` (335086)
<pre>
BAM_INPUT
TOTAL READS     = 13139331
COUNTS_SUM      = 13139331
DISTINCT READS  = 9.95264e+06
DISTINCT COUNTS = 149
MAX COUNT       = 430
COUNTS OF 1     = 8.33143e+06
OBSERVED COUNTS (431)
1	8331431
2	1064073
3	273280
4	110730
5	57100
6	33562
7	21886
8	14351
9	10164
10	7277
11	5449
12	4060
13	3259
14	2572
15	2117
16	1690
17	1348
18	1122
19	924
20	809
21	697
22	515
23	473
24	408
25	362
26	314
27	268
28	197
29	186
30	150
31	147
32	163
33	117
34	133
35	90
36	92
37	78
38	88
39	59
40	50
41	60
42	49
43	34
44	46
45	41
46	32
47	44
48	36
49	27
50	22
51	33
52	22
53	21
54	16
55	20
56	17
57	20
58	12
59	14
60	10
61	9
62	8
63	12
64	9
65	11
66	6
67	8
68	5
69	10
70	8
71	8
72	7
73	3
74	8
75	3
76	6
77	5
78	2
79	3
80	6
81	5
82	6
83	7
84	3
85	4
86	2
87	2
88	3
89	3
90	4
91	2
92	6
93	2
94	4
95	3
96	2
98	1
99	2
100	2
101	1
102	3
104	1
105	1
107	2
108	2
109	1
110	5
111	2
112	1
113	1
114	1
116	2
117	2
118	2
123	1
125	1
126	1
128	1
132	1
133	2
134	1
135	1
136	1
139	1
141	1
143	1
144	1
147	2
150	1
151	1
156	2
167	2
168	1
172	2
174	1
177	1
189	1
191	1
193	2
215	1
216	1
225	1
280	1
296	1
349	1
361	1
388	1
402	1
430	1

sample size: 1000000
sample size: 2000000
sample size: 3000000
sample size: 4000000
sample size: 5000000
sample size: 6000000
sample size: 7000000
sample size: 8000000
sample size: 9000000
sample size: 10000000
sample size: 11000000
sample size: 12000000
sample size: 13000000
</pre>
Command completed. Elapsed time: 0:01:14. Running peak memory: 3.606GB.  
  PID: 335086;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort_dups.bam` (335365)
<pre>
BAM_INPUT
TOTAL READS     = 13139331
DISTINCT READS  = 9.95264e+06
DISTINCT COUNTS = 149
MAX COUNT       = 430
COUNTS OF 1     = 8.33143e+06
MAX TERMS       = 96
OBSERVED COUNTS (431)
1	8331431
2	1064073
3	273280
4	110730
5	57100
6	33562
7	21886
8	14351
9	10164
10	7277
11	5449
12	4060
13	3259
14	2572
15	2117
16	1690
17	1348
18	1122
19	924
20	809
21	697
22	515
23	473
24	408
25	362
26	314
27	268
28	197
29	186
30	150
31	147
32	163
33	117
34	133
35	90
36	92
37	78
38	88
39	59
40	50
41	60
42	49
43	34
44	46
45	41
46	32
47	44
48	36
49	27
50	22
51	33
52	22
53	21
54	16
55	20
56	17
57	20
58	12
59	14
60	10
61	9
62	8
63	12
64	9
65	11
66	6
67	8
68	5
69	10
70	8
71	8
72	7
73	3
74	8
75	3
76	6
77	5
78	2
79	3
80	6
81	5
82	6
83	7
84	3
85	4
86	2
87	2
88	3
89	3
90	4
91	2
92	6
93	2
94	4
95	3
96	2
98	1
99	2
100	2
101	1
102	3
104	1
105	1
107	2
108	2
109	1
110	5
111	2
112	1
113	1
114	1
116	2
117	2
118	2
123	1
125	1
126	1
128	1
132	1
133	2
134	1
135	1
136	1
139	1
141	1
143	1
144	1
147	2
150	1
151	1
156	2
167	2
168	1
172	2
174	1
177	1
189	1
191	1
193	2
215	1
216	1
225	1
280	1
296	1
349	1
361	1
388	1
402	1
430	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
................_.._.............................................._....................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:14. Running peak memory: 3.606GB.  
  PID: 335365;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort_dups.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam) > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_counts.txt` (335432)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 3.606GB.  
  PID: 335432;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_plot` (335456)
<pre>
Processing Jurkat_ChRO-seq_1
INFO: Found real counts for Jurkat_ChRO-seq_1 - Total (M): 21.334642 Unique (M): 13.139331

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.606GB.  
  PID: 335456;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Library complexity`	QC_hg38/Jurkat_ChRO-seq_1_preseq_plot.pdf	Library complexity	QC_hg38/Jurkat_ChRO-seq_1_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.7897	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-18 12:07:45) elapsed: 177.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -c 8 -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_bamQC.tsv` (335477)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/tmp_Jurkat_ChRO-seq_1_sort_2cs7e6bp'
Processing with 8 cores...
Discarding 109 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr22_KI270732v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1']
Keeping 86 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.606GB.  
  PID: 335477;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.229GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_bamQC.tsv`

> `NRF`	0.38	PEPPRO	_RES_

> `PBC1`	0.64	PEPPRO	_RES_

> `PBC2`	3.39	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_unmap.bam`  

> `samtools view -b -@ 8 -f 4  /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_unmap.bam` (335559)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.606GB.  
  PID: 335559;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp.bam`

> `Unmapped_reads`	1113286	PEPPRO	_RES_

### Split BAM by strand (02-18 12:08:22) elapsed: 37.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam` (335589)
<pre>
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 3.606GB.  
  PID: 335589;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam` (335650)
<pre>
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 3.606GB.  
  PID: 335650;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-18 12:10:24) elapsed: 122.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (335958)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.606GB.  
  PID: 335958;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_plus_TssEnrichment.txt` (335959)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.606GB.  
  PID: 335959;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.371GB


> `TSS_Plus_Score`	55.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_minus_TssEnrichment.txt` (335987)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.606GB.  
  PID: 335987;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.378GB


> `TSS_Minus_Score`	12.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_minus_TssEnrichment.txt` (336014)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.606GB.  
  PID: 336014;	Command: Rscript;	Return code: 0;	Memory used: 0.235GB

> `TSS enrichment`	QC_hg38/Jurkat_ChRO-seq_1_TSSenrichment.pdf	TSS enrichment	QC_hg38/Jurkat_ChRO-seq_1_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt` (336033,336034,336035,336036)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.606GB.  
  PID: 336033;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 336035;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 336034;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 336036;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_keep.txt` (336038)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.606GB.  
  PID: 336038;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (02-18 12:10:50) elapsed: 26.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_tss.bed` (336040,336041)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.606GB.  
  PID: 336040;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 336041;	Command: bedtools;	Return code: 0;	Memory used: 0.096GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed` (336045,336046)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.606GB.  
  PID: 336045;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 336046;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_TSS_density.bed` (336048,336049,336050,336051)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 3.606GB.  
  PID: 336049;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 336051;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 336048;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB  
  PID: 336050;	Command: sort;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_gene_body_density.bed` (336081,336082,336083)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 3.606GB.  
  PID: 336083;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 336081;	Command: bedtools;	Return code: 0;	Memory used: 0.098GB  
  PID: 336082;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed` (336115,336116,336117)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.606GB.  
  PID: 336115;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 336117;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 336116;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	50.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed` (336123)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.606GB.  
  PID: 336123;	Command: Rscript;	Return code: 0;	Memory used: 0.221GB

> `Pause index`	QC_hg38/Jurkat_ChRO-seq_1_pause_index.pdf	Pause index	QC_hg38/Jurkat_ChRO-seq_1_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed` (336145)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.606GB.  
  PID: 336145;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Calculate Fraction of Reads in pre-mature mRNA (02-18 12:12:05) elapsed: 75.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam`
21338079 7688297

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam`
21338079 7358614

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_gene_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_gene_sort.bed` (336210,336211)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.606GB.  
  PID: 336210;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 336211;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_gene_coverage.bed` (336213)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 3.606GB.  
  PID: 336213;	Command: bedtools;	Return code: 0;	Memory used: 0.075GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/raw/hg38_annotations.bed.gz` (336253)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.606GB.  
  PID: 336253;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/raw/hg38_annotations.bed` (336254)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.606GB.  
  PID: 336254;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate fraction and proportion of reads in features (FRiF/PRiF) (02-18 12:13:28) elapsed: 83.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/raw/hg38_annotations.bed` (336263)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.606GB.  
  PID: 336263;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3_UTR"` (336265)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.606GB.  
  PID: 336265;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3_UTR_sort.bed` (336266,336267,336268,336269)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.606GB.  
  PID: 336266;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 336267;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 336269;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 336268;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_3_UTR_plus_coverage.bed` (336271)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.606GB.  
  PID: 336271;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_3_UTR_minus_coverage.bed` (336284)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.606GB.  
  PID: 336284;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5_UTR"` (336296)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.606GB.  
  PID: 336296;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5_UTR_sort.bed` (336297,336298,336299,336300)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.606GB.  
  PID: 336297;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 336298;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 336300;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 336299;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_5_UTR_plus_coverage.bed` (336302)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.606GB.  
  PID: 336302;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_5_UTR_minus_coverage.bed` (336319)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.606GB.  
  PID: 336319;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Enhancer_sort.bed` (336331,336332,336333,336334)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.606GB.  
  PID: 336331;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 336332;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 336334;	Command: bedtools;	Return code: 0;	Memory used: 0.047GB  
  PID: 336333;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Enhancer_plus_coverage.bed` (336337)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.606GB.  
  PID: 336337;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Enhancer_minus_coverage.bed` (336349)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.606GB.  
  PID: 336349;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Exon_sort.bed` (336361,336362,336363,336364)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.606GB.  
  PID: 336361;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 336362;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 336364;	Command: bedtools;	Return code: 0;	Memory used: 0.176GB  
  PID: 336363;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Exon_plus_coverage.bed` (336369)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.606GB.  
  PID: 336369;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Exon_minus_coverage.bed` (336599)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.606GB.  
  PID: 336599;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Intron_sort.bed` (336622,336623,336624,336625)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.606GB.  
  PID: 336622;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 336624;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 336623;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 336625;	Command: bedtools;	Return code: 0;	Memory used: 0.083GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Intron_plus_coverage.bed` (336628)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.606GB.  
  PID: 336628;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Intron_minus_coverage.bed` (336643)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.606GB.  
  PID: 336643;	Command: bedtools;	Return code: 0;	Memory used: 0.047GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_sort.bed` (336658,336659,336660,336661)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.606GB.  
  PID: 336658;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 336659;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 336661;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 336660;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_plus_coverage.bed` (336663)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.606GB.  
  PID: 336663;	Command: bedtools;	Return code: 0;	Memory used: 0.055GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_minus_coverage.bed` (336681)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.606GB.  
  PID: 336681;	Command: bedtools;	Return code: 0;	Memory used: 0.038GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_Flanking_Region"` (336697)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.606GB.  
  PID: 336697;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed` (336698,336699,336700,336701)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.606GB.  
  PID: 336698;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 336700;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 336699;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 336701;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (336704)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.606GB.  
  PID: 336704;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed` (336718)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.606GB.  
  PID: 336718;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB


### Plot FRiF/PRiF (02-18 12:16:55) elapsed: 207.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s Jurkat_ChRO-seq_1 -z 3099922541 -n 10780321 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_frif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (336741)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 3.606GB.  
  PID: 336741;	Command: Rscript;	Return code: 0;	Memory used: 0.448GB

> `FRiF`	QC_hg38/Jurkat_ChRO-seq_1_frif.pdf	FRiF	QC_hg38/Jurkat_ChRO-seq_1_frif.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_prif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s Jurkat_ChRO-seq_1 -z 3099922541 -n 10780321 -y prif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_prif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (336787)
<pre>
Cumulative prif plot completed!

</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.606GB.  
  PID: 336787;	Command: Rscript;	Return code: 0;	Memory used: 0.466GB

> `PRiF`	QC_hg38/Jurkat_ChRO-seq_1_prif.pdf	PRiF	QC_hg38/Jurkat_ChRO-seq_1_prif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-18 12:18:01) elapsed: 67.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_exons_sort.bed` (336821,336822)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.606GB.  
  PID: 336822;	Command: bedtools;	Return code: 0;	Memory used: 0.09GB  
  PID: 336821;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_introns_sort.bed` (336831,336832,336833)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.606GB.  
  PID: 336831;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 336833;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 336832;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exons_coverage.bed` (336840)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.606GB.  
  PID: 336840;	Command: bedtools;	Return code: 0;	Memory used: 0.023GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_introns_coverage.bed` (336864)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.606GB.  
  PID: 336864;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/21.338079)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exons_rpkm.bed` (336897,336898,336899)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.606GB.  
  PID: 336897;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 336899;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 336898;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/21.338079)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_introns_rpkm.bed` (336901,336902,336903)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.606GB.  
  PID: 336901;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 336903;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 336902;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exon_intron_ratios.bed` (336906,336907,336908)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.606GB.  
  PID: 336906;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 336908;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 336907;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.19	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exon_intron_ratios.bed --annotate` (336914)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.606GB.  
  PID: 336914;	Command: Rscript;	Return code: 0;	Memory used: 0.239GB

> `mRNA contamination`	QC_hg38/Jurkat_ChRO-seq_1_mRNA_contamination.pdf	mRNA contamination	QC_hg38/Jurkat_ChRO-seq_1_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exon_intron_ratios.bed` (336933)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.606GB.  
  PID: 336933;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (02-18 12:19:19) elapsed: 77.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam` (336941)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.606GB.  
  PID: 336941;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (336948)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam'
Temporary files will be stored in: 'tmp_Jurkat_ChRO-seq_1_plus_cuttrace_6p4y0zuw'
Processing with 2 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 122 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr9_KI270717v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr15_KI270727v1_random', 'chr22_KI270732v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270743v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270742v1']
Keeping 73 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270579v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_GL000213v1', 'chrUn_KI270744v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 73 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_plus_exact_body_0-mer.bw'
Merging 73 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:42. Running peak memory: 3.606GB.  
  PID: 336948;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.671GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam` (338627)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.606GB.  
  PID: 338627;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_minus_exact_body_0-mer.bw -p 5 --variable-step --tail-edge` (338634)
<pre>
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam'
Temporary files will be stored in: 'tmp_Jurkat_ChRO-seq_1_minus_cuttrace_4dlcm764'
Processing with 5 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 126 chunk(s) of reads: ['chrM', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr14_KI270725v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1']
Keeping 69 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270538v1', 'chrUn_KI270580v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 69 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_minus_exact_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.606GB.  
  PID: 338634;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.262GB

Starting cleanup: 63 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  1:20:37
*  Total elapsed time (all runs):  1:56:41
*         Peak memory (this run):  3.606 GB
*        Pipeline completed time: 2020-02-18 12:29:34
