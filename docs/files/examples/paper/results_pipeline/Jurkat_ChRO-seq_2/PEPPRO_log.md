### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name Jurkat_ChRO-seq_2 --genome hg38 --input /project/shefflab/data/sra_fastq/SRR7616134.fastq.gz --single-or-paired single --protocol PRO --umi-len 6 --max-len -1 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-aj38-13c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/
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
*              `input`:  `['/project/shefflab/data/sra_fastq/SRR7616134.fastq.gz']`
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
*        `sample_name`:  `Jurkat_ChRO-seq_2`
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

Local input file: /project/shefflab/data/sra_fastq/SRR7616134.fastq.gz

> `File_mb`	1914.91	PEPPRO	_RES_

> `Read_type`	single	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-18 11:08:58) elapsed: 0.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/raw/Jurkat_ChRO-seq_2.fastq.gz`  

> `ln -sf /project/shefflab/data/sra_fastq/SRR7616134.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/raw/Jurkat_ChRO-seq_2.fastq.gz` (327687)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 327687;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/raw/Jurkat_ChRO-seq_2.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/raw/Jurkat_ChRO-seq_2.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1.fastq` (327688)
<pre>
</pre>
Command completed. Elapsed time: 0:02:41. Running peak memory: 0.003GB.  
  PID: 327688;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	49841170	PEPPRO	_RES_

> `Fastq_reads`	49841170	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/raw/Jurkat_ChRO-seq_2.fastq.gz']

### FASTQ processing:  (02-18 11:12:10) elapsed: 192.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1_processed.fastq`  

> `(cutadapt -j 8 -m 8 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/cutadapt/Jurkat_ChRO-seq_2_R1_cutadapt.txt` (328152)
<pre>
</pre>
Command completed. Elapsed time: 0:01:50. Running peak memory: 0.293GB.  
  PID: 328152;	Command: cutadapt;	Return code: 0;	Memory used: 0.293GB


> `seqtk trimfq -b 6 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1_noadap.fastq | seqtk seq -L 8 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1_processed.fastq` (328306,328307)
<pre>
</pre>
Command completed. Elapsed time: 0:00:59. Running peak memory: 0.293GB.  
  PID: 328306;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 328307;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	48264585	PEPPRO	_RES_

> `Trim_loss_rate`	3.16	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1_processed.fastq` (329739)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 5% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 10% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 15% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 20% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 25% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 30% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 35% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 40% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 45% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 50% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 55% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 60% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 65% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 70% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 75% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 80% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 85% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 90% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Approx 95% complete for Jurkat_ChRO-seq_2_R1_processed.fastq
Analysis complete for Jurkat_ChRO-seq_2_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:02:25. Running peak memory: 0.293GB.  
  PID: 329739;	Command: fastqc;	Return code: 0;	Memory used: 0.184GB

> `FastQC report r1`	fastqc/Jurkat_ChRO-seq_2_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1_dedup.fastq /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1_noadap.fastq` (329992)
<pre>
[INFO][0m 5764233 duplicated records removed
</pre>
Command completed. Elapsed time: 0:01:44. Running peak memory: 3.993GB.  
  PID: 329992;	Command: seqkit;	Return code: 0;	Memory used: 3.993GB


> `seqtk trimfq -b 6 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1_dedup.fastq | seqtk seq -L 8 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1_trimmed.fastq` (330096,330097)
<pre>
</pre>
Command completed. Elapsed time: 0:00:56. Running peak memory: 3.993GB.  
  PID: 330096;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 330097;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/cutadapt/Jurkat_ChRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	42416468.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/cutadapt/Jurkat_ChRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/cutadapt/Jurkat_ChRO-seq_2_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/cutadapt/Jurkat_ChRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_too_short`	223715.0	PEPPRO	_RES_

> `Duplicate_reads`	5764233.0	PEPPRO	_RES_

> `Pct_reads_too_short`	0.4489	PEPPRO	_RES_

### Plot adapter insertion distribution (02-18 11:20:43) elapsed: 513.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/cutadapt/Jurkat_ChRO-seq_2_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/cutadapt/Jurkat_ChRO-seq_2_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/cutadapt -u 6` (330433)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.993GB.  
  PID: 330433;	Command: Rscript;	Return code: 0;	Memory used: 0.125GB

> `Adapter insertion distribution`	cutadapt/Jurkat_ChRO-seq_2_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/Jurkat_ChRO-seq_2_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/cutadapt/Jurkat_ChRO-seq_2_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	39	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-18 11:20:51) elapsed: 8.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/cutadapt/Jurkat_ChRO-seq_2_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/cutadapt/Jurkat_ChRO-seq_2_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/cutadapt/Jurkat_ChRO-seq_2_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/cutadapt/Jurkat_ChRO-seq_2_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/cutadapt/Jurkat_ChRO-seq_2_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.4899	PEPPRO	_RES_

### Prealignments (02-18 11:20:51) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-18 11:20:51) elapsed: 0.0 _TIME_


> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id Jurkat_ChRO-seq_2 -U /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/prealignments/Jurkat_ChRO-seq_2_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
48264585 reads; of these:
  48264585 (100.00%) were unpaired; of these:
    44313938 (91.81%) aligned 0 times
    3950647 (8.19%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
8.19% overall alignment rate

> `Aligned_reads_human_rDNA`	3950647.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	8.19	PEPPRO	_RES_

### Map to human_rDNA (02-18 11:25:44) elapsed: 293.0 _TIME_


> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id Jurkat_ChRO-seq_2 -U /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/fastq/Jurkat_ChRO-seq_2_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/prealignments/Jurkat_ChRO-seq_2_human_rDNA_unmap_dups.fq 2>&1 > /dev/null)`

### Map to genome (02-18 11:29:34) elapsed: 231.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id Jurkat_ChRO-seq_2 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/prealignments/Jurkat_ChRO-seq_2_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/tmp3yfbs809 -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_temp.bam` (331215,331216,331217)
<pre>
44313938 reads; of these:
  44313938 (100.00%) were unpaired; of these:
    2435871 (5.50%) aligned 0 times
    29426081 (66.40%) aligned exactly 1 time
    12451986 (28.10%) aligned >1 times
94.50% overall alignment rate
[bam_sort_core] merging from 11 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:26:35. Running peak memory: 3.993GB.  
  PID: 331215;	Command: bowtie2;	Return code: 0;	Memory used: 3.637GB  
  PID: 331216;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 331217;	Command: samtools;	Return code: 0;	Memory used: 0.874GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam` (334291)
<pre>
</pre>
Command completed. Elapsed time: 0:01:14. Running peak memory: 3.993GB.  
  PID: 334291;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	41878067	PEPPRO	_RES_

> `QC_filtered_reads`	7874013	PEPPRO	_RES_

> `Aligned_reads`	34004054	PEPPRO	_RES_

> `Alignment_rate`	70.45	PEPPRO	_RES_

> `Total_efficiency`	68.22	PEPPRO	_RES_

> `Read_depth`	4.27	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id Jurkat_ChRO-seq_2 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/prealignments/Jurkat_ChRO-seq_2_human_rDNA_unmap_dups.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/tmp3yfbs809 -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_temp_dups.bam` (335678,335679,335681)
<pre>
39606611 reads; of these:
  39606611 (100.00%) were unpaired; of these:
    2278810 (5.75%) aligned 0 times
    26411590 (66.68%) aligned exactly 1 time
    10916211 (27.56%) aligned >1 times
94.25% overall alignment rate
[bam_sort_core] merging from 10 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:22:11. Running peak memory: 3.993GB.  
  PID: 335678;	Command: bowtie2;	Return code: 0;	Memory used: 3.61GB  
  PID: 335679;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 335681;	Command: samtools;	Return code: 0;	Memory used: 0.931GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_temp_dups.bam > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort_dups.bam` (339372)
<pre>
</pre>
Command completed. Elapsed time: 0:01:10. Running peak memory: 3.993GB.  
  PID: 339372;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


### Compress all unmapped read files (02-18 12:33:16) elapsed: 3822.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/prealignments/Jurkat_ChRO-seq_2_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/prealignments/Jurkat_ChRO-seq_2_human_rDNA_unmap.fq` (339457)
<pre>
</pre>
Command completed. Elapsed time: 0:01:14. Running peak memory: 3.993GB.  
  PID: 339457;	Command: pigz;	Return code: 0;	Memory used: 0.011GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_temp.bam` (339531)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 3.993GB.  
  PID: 339531;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	18250	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam` (339560)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 3.993GB.  
  PID: 339560;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_noMT.bam` (339797,339798,339799,339800)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 3.993GB.  
  PID: 339798;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 339797;	Command: samtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 339799;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 339800;	Command: xargs;	Return code: 0;	Memory used: 0.064GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam` (339838)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.993GB.  
  PID: 339838;	Command: mv;	Return code: 0;	Memory used: 0.0GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam` (339840)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 3.993GB.  
  PID: 339840;	Command: samtools;	Return code: 0;	Memory used: 0.012GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	70	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_temp_dups.bam` (339926)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.993GB.  
  PID: 339926;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort_dups.bam`  

### Calculate library complexity (02-18 12:37:51) elapsed: 275.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_preseq_out.txt -B /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort_dups.bam` (339953)
<pre>
BAM_INPUT
TOTAL READS     = 30468400
COUNTS_SUM      = 30468400
DISTINCT READS  = 1.6294e+07
DISTINCT COUNTS = 1118
MAX COUNT       = 6926
COUNTS OF 1     = 1.31032e+07
OBSERVED COUNTS (6927)
1	13103248
2	1726484
3	540128
4	254315
5	147467
6	96195
7	68125
8	50567
9	38638
10	30601
11	24915
12	20508
13	17467
14	14643
15	12769
16	11038
17	9701
18	8575
19	7601
20	6875
21	6210
22	5518
23	5088
24	4592
25	4188
26	3754
27	3445
28	3289
29	3024
30	2717
31	2588
32	2388
33	2274
34	2083
35	1899
36	1954
37	1666
38	1757
39	1529
40	1529
41	1388
42	1291
43	1263
44	1190
45	1148
46	1044
47	1022
48	1000
49	892
50	889
51	866
52	835
53	786
54	762
55	698
56	720
57	663
58	663
59	576
60	563
61	564
62	562
63	446
64	491
65	459
66	481
67	422
68	461
69	374
70	387
71	407
72	394
73	365
74	388
75	349
76	317
77	319
78	327
79	317
80	292
81	306
82	282
83	283
84	271
85	249
86	243
87	263
88	248
89	234
90	226
91	209
92	214
93	208
94	191
95	187
96	197
97	195
98	194
99	182
100	187
101	170
102	173
103	179
104	162
105	134
106	142
107	138
108	143
109	124
110	141
111	156
112	146
113	122
114	126
115	121
116	124
117	113
118	117
119	121
120	107
121	116
122	129
123	104
124	107
125	103
126	94
127	90
128	79
129	94
130	103
131	94
132	97
133	92
134	78
135	84
136	98
137	73
138	78
139	73
140	76
141	82
142	77
143	62
144	71
145	84
146	65
147	79
148	62
149	69
150	73
151	63
152	65
153	64
154	59
155	62
156	65
157	56
158	60
159	46
160	60
161	50
162	47
163	55
164	47
165	41
166	54
167	55
168	49
169	49
170	48
171	50
172	64
173	52
174	36
175	56
176	42
177	45
178	42
179	40
180	37
181	32
182	57
183	46
184	41
185	39
186	43
187	37
188	37
189	35
190	32
191	37
192	48
193	31
194	31
195	36
196	34
197	32
198	32
199	32
200	37
201	35
202	31
203	30
204	22
205	38
206	32
207	34
208	27
209	15
210	22
211	18
212	29
213	32
214	30
215	33
216	32
217	36
218	17
219	20
220	23
221	19
222	36
223	22
224	26
225	28
226	25
227	26
228	28
229	18
230	23
231	29
232	24
233	20
234	28
235	30
236	23
237	14
238	19
239	16
240	10
241	9
242	22
243	29
244	17
245	18
246	11
247	17
248	23
249	15
250	19
251	13
252	16
253	28
254	20
255	14
256	16
257	11
258	14
259	13
260	9
261	16
262	15
263	14
264	16
265	15
266	12
267	17
268	17
269	20
270	11
271	15
272	17
273	12
274	15
275	15
276	29
277	15
278	14
279	22
280	11
281	13
282	19
283	9
284	15
285	17
286	9
287	16
288	16
289	10
290	8
291	15
292	10
293	8
294	10
295	14
296	11
297	12
298	14
299	3
300	14
301	12
302	12
303	14
304	5
305	18
306	10
307	9
308	8
309	16
310	13
311	16
312	10
313	11
314	14
315	11
316	16
317	7
318	9
319	10
320	10
321	7
322	6
323	9
324	10
325	11
326	10
327	12
328	16
329	6
330	8
331	6
332	10
333	7
334	9
335	7
336	10
337	5
338	8
339	5
340	5
341	5
342	6
343	6
344	9
345	3
346	13
347	11
348	2
349	5
350	10
351	7
352	6
353	6
354	6
355	6
356	8
357	9
358	8
359	9
360	2
361	6
362	7
363	6
364	2
365	2
366	5
367	6
368	6
369	7
370	4
371	3
372	3
373	10
374	8
375	7
376	5
377	4
378	5
379	3
380	6
381	3
382	4
383	6
384	4
385	1
386	5
387	3
388	5
389	8
390	3
391	5
392	3
393	5
394	2
395	7
396	3
397	6
398	7
399	3
400	3
401	3
402	5
403	3
404	7
405	6
406	1
407	7
408	7
409	6
410	5
411	7
412	2
413	6
414	9
415	3
416	3
417	4
418	5
419	3
420	6
421	1
422	4
423	5
424	7
425	4
426	5
427	7
428	4
429	7
430	6
431	11
432	7
433	3
434	5
435	7
436	4
437	3
438	6
439	4
440	2
441	5
442	6
443	5
444	3
445	3
446	3
447	5
448	4
449	7
450	2
451	5
452	3
453	6
454	2
455	4
456	5
457	8
458	4
459	5
460	3
461	3
462	1
463	4
464	3
465	4
466	1
467	2
468	6
469	4
470	4
471	3
472	3
473	4
474	2
475	2
476	4
477	1
478	6
479	4
480	1
481	5
482	1
483	4
484	2
486	4
487	2
488	5
489	4
490	6
491	3
492	5
493	8
494	2
495	3
496	1
497	2
498	5
499	1
500	4
501	6
502	4
503	1
504	1
505	4
506	4
507	1
508	4
509	2
510	3
511	2
513	6
514	2
515	2
516	6
517	7
518	2
519	3
520	1
521	6
522	2
523	5
524	1
525	2
526	2
527	2
528	1
529	3
530	2
531	2
532	1
533	4
534	8
535	1
536	2
537	4
538	2
539	2
540	1
541	1
542	1
543	1
544	2
545	3
546	1
547	2
548	4
549	2
550	1
552	2
553	3
554	3
555	3
556	5
557	4
558	2
559	4
561	3
562	4
563	1
564	4
567	3
569	2
570	2
571	1
572	3
573	1
574	6
575	3
576	2
577	2
578	2
579	1
580	4
581	4
582	3
583	1
584	3
585	3
586	1
587	2
588	2
589	2
590	1
591	3
592	1
593	2
594	6
595	4
596	6
597	3
598	2
599	2
601	5
603	2
604	5
605	2
606	2
607	3
608	2
609	5
610	1
611	2
612	2
613	4
616	1
617	2
618	1
619	1
620	2
621	2
622	1
623	1
624	3
625	5
626	1
627	2
628	4
630	1
631	2
632	2
633	4
635	4
636	3
637	2
638	2
639	4
640	3
641	1
642	2
644	1
645	4
646	4
647	2
648	1
649	2
650	2
651	2
652	3
653	2
654	4
655	1
656	2
657	1
658	2
659	3
660	1
661	5
662	1
663	3
664	3
665	1
666	2
667	1
668	1
671	1
672	2
673	1
674	1
675	1
677	2
678	2
679	3
680	1
682	4
684	2
685	1
686	1
687	1
688	2
689	1
690	1
691	1
694	4
695	2
696	2
697	1
699	2
700	2
702	1
703	1
704	1
705	1
706	1
707	1
712	1
713	3
714	1
715	1
716	2
717	1
718	1
719	2
720	2
721	3
722	2
723	1
724	2
725	1
729	1
730	1
731	1
732	3
733	2
735	1
736	2
737	1
738	1
739	1
741	1
744	3
746	2
747	1
748	1
749	1
750	1
751	1
752	1
753	2
754	1
755	1
756	1
760	1
761	1
763	3
764	1
765	1
766	1
767	1
768	2
770	3
771	4
772	1
774	1
775	1
777	1
778	1
779	1
780	2
782	1
784	2
789	1
790	1
792	2
793	4
794	1
797	1
800	1
802	2
803	1
804	1
806	2
807	1
808	3
809	2
810	3
811	3
813	1
814	1
815	1
816	1
817	1
818	2
819	1
820	1
821	1
822	1
823	1
824	3
828	1
829	2
832	1
835	3
836	1
837	2
838	1
840	1
844	1
847	2
849	1
850	1
852	1
854	3
858	1
859	2
861	1
862	1
863	3
864	1
866	1
867	2
869	2
871	1
872	1
873	2
876	4
878	2
880	1
882	2
885	2
886	1
887	1
889	1
890	1
891	1
892	1
895	2
897	1
898	2
899	1
902	2
903	1
905	2
907	2
909	2
910	1
911	1
912	1
917	1
918	1
921	1
922	2
923	2
927	2
928	1
930	1
932	1
937	1
939	1
940	1
944	1
945	1
947	1
951	2
952	1
953	2
956	2
959	1
961	1
963	1
964	2
965	2
967	1
970	1
971	1
981	2
982	2
984	2
985	1
986	1
989	1
990	1
992	1
993	2
994	1
995	1
996	1
997	1
998	1
1000	2
1002	1
1003	1
1004	2
1005	1
1012	1
1013	1
1017	1
1021	2
1022	1
1025	1
1027	1
1030	1
1031	2
1033	2
1034	1
1036	1
1037	1
1038	2
1040	1
1041	1
1046	1
1051	1
1053	1
1056	1
1058	1
1061	1
1062	1
1063	2
1065	1
1067	1
1069	2
1070	1
1072	1
1073	1
1074	1
1075	1
1077	1
1079	1
1082	1
1083	1
1084	2
1086	1
1093	1
1095	2
1096	1
1099	1
1101	1
1109	1
1110	1
1115	1
1122	1
1124	2
1125	1
1127	1
1130	1
1132	1
1139	2
1140	2
1142	1
1143	1
1149	1
1151	1
1153	2
1156	1
1158	1
1162	1
1163	1
1164	1
1167	2
1168	1
1170	1
1177	1
1180	1
1183	1
1185	1
1189	2
1200	1
1202	1
1203	1
1204	1
1206	1
1212	2
1215	1
1216	1
1217	2
1218	1
1219	1
1224	1
1225	1
1228	1
1231	1
1238	1
1243	1
1246	1
1249	1
1259	1
1264	1
1269	2
1274	1
1283	1
1284	1
1286	1
1288	1
1289	1
1295	1
1305	1
1313	1
1317	1
1319	1
1326	1
1328	1
1341	1
1345	1
1350	1
1351	1
1354	2
1359	1
1365	1
1371	1
1375	1
1386	1
1388	1
1394	1
1422	1
1426	1
1427	1
1433	1
1435	1
1438	1
1439	1
1445	1
1452	1
1456	2
1462	1
1465	1
1472	1
1476	1
1478	1
1485	1
1486	1
1487	1
1492	1
1493	1
1499	1
1516	1
1533	1
1535	1
1536	1
1549	1
1550	1
1556	1
1558	1
1569	1
1588	2
1609	1
1629	1
1630	1
1644	2
1649	1
1652	1
1660	1
1686	2
1689	1
1692	1
1703	1
1712	1
1713	1
1719	1
1727	1
1730	1
1738	1
1741	1
1757	1
1759	1
1762	1
1772	1
1792	1
1824	1
1831	1
1833	1
1847	1
1855	1
1870	1
1878	1
1880	1
1900	1
1920	1
1924	1
1927	1
1934	2
1941	1
1970	1
1977	1
1987	1
1991	1
1999	1
2013	1
2016	1
2034	1
2036	1
2046	1
2047	1
2050	1
2054	1
2072	1
2076	1
2103	1
2108	1
2152	1
2154	1
2173	1
2176	1
2184	1
2186	1
2188	1
2236	1
2273	1
2282	1
2296	1
2414	1
2417	1
2443	1
2489	1
2548	1
2558	1
2570	2
2586	1
2590	1
2676	1
2690	1
2697	2
2714	1
2732	1
2736	1
2859	1
2900	1
2906	1
2935	1
2941	1
2965	1
2969	1
3049	1
3104	1
3109	1
3177	1
3432	1
3478	1
3536	1
3644	1
3663	1
4052	1
4500	1
4782	1
5341	1
5798	1
5871	1
5948	1
6559	1
6926	1

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
sample size: 14000000
sample size: 15000000
sample size: 16000000
sample size: 17000000
sample size: 18000000
sample size: 19000000
sample size: 20000000
sample size: 21000000
sample size: 22000000
sample size: 23000000
sample size: 24000000
sample size: 25000000
sample size: 26000000
sample size: 27000000
sample size: 28000000
sample size: 29000000
sample size: 30000000
</pre>
Command completed. Elapsed time: 0:02:47. Running peak memory: 3.993GB.  
  PID: 339953;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort_dups.bam` (340488)
<pre>
BAM_INPUT
TOTAL READS     = 30468400
DISTINCT READS  = 1.6294e+07
DISTINCT COUNTS = 1118
MAX COUNT       = 6926
COUNTS OF 1     = 1.31032e+07
MAX TERMS       = 100
OBSERVED COUNTS (6927)
1	13103248
2	1726484
3	540128
4	254315
5	147467
6	96195
7	68125
8	50567
9	38638
10	30601
11	24915
12	20508
13	17467
14	14643
15	12769
16	11038
17	9701
18	8575
19	7601
20	6875
21	6210
22	5518
23	5088
24	4592
25	4188
26	3754
27	3445
28	3289
29	3024
30	2717
31	2588
32	2388
33	2274
34	2083
35	1899
36	1954
37	1666
38	1757
39	1529
40	1529
41	1388
42	1291
43	1263
44	1190
45	1148
46	1044
47	1022
48	1000
49	892
50	889
51	866
52	835
53	786
54	762
55	698
56	720
57	663
58	663
59	576
60	563
61	564
62	562
63	446
64	491
65	459
66	481
67	422
68	461
69	374
70	387
71	407
72	394
73	365
74	388
75	349
76	317
77	319
78	327
79	317
80	292
81	306
82	282
83	283
84	271
85	249
86	243
87	263
88	248
89	234
90	226
91	209
92	214
93	208
94	191
95	187
96	197
97	195
98	194
99	182
100	187
101	170
102	173
103	179
104	162
105	134
106	142
107	138
108	143
109	124
110	141
111	156
112	146
113	122
114	126
115	121
116	124
117	113
118	117
119	121
120	107
121	116
122	129
123	104
124	107
125	103
126	94
127	90
128	79
129	94
130	103
131	94
132	97
133	92
134	78
135	84
136	98
137	73
138	78
139	73
140	76
141	82
142	77
143	62
144	71
145	84
146	65
147	79
148	62
149	69
150	73
151	63
152	65
153	64
154	59
155	62
156	65
157	56
158	60
159	46
160	60
161	50
162	47
163	55
164	47
165	41
166	54
167	55
168	49
169	49
170	48
171	50
172	64
173	52
174	36
175	56
176	42
177	45
178	42
179	40
180	37
181	32
182	57
183	46
184	41
185	39
186	43
187	37
188	37
189	35
190	32
191	37
192	48
193	31
194	31
195	36
196	34
197	32
198	32
199	32
200	37
201	35
202	31
203	30
204	22
205	38
206	32
207	34
208	27
209	15
210	22
211	18
212	29
213	32
214	30
215	33
216	32
217	36
218	17
219	20
220	23
221	19
222	36
223	22
224	26
225	28
226	25
227	26
228	28
229	18
230	23
231	29
232	24
233	20
234	28
235	30
236	23
237	14
238	19
239	16
240	10
241	9
242	22
243	29
244	17
245	18
246	11
247	17
248	23
249	15
250	19
251	13
252	16
253	28
254	20
255	14
256	16
257	11
258	14
259	13
260	9
261	16
262	15
263	14
264	16
265	15
266	12
267	17
268	17
269	20
270	11
271	15
272	17
273	12
274	15
275	15
276	29
277	15
278	14
279	22
280	11
281	13
282	19
283	9
284	15
285	17
286	9
287	16
288	16
289	10
290	8
291	15
292	10
293	8
294	10
295	14
296	11
297	12
298	14
299	3
300	14
301	12
302	12
303	14
304	5
305	18
306	10
307	9
308	8
309	16
310	13
311	16
312	10
313	11
314	14
315	11
316	16
317	7
318	9
319	10
320	10
321	7
322	6
323	9
324	10
325	11
326	10
327	12
328	16
329	6
330	8
331	6
332	10
333	7
334	9
335	7
336	10
337	5
338	8
339	5
340	5
341	5
342	6
343	6
344	9
345	3
346	13
347	11
348	2
349	5
350	10
351	7
352	6
353	6
354	6
355	6
356	8
357	9
358	8
359	9
360	2
361	6
362	7
363	6
364	2
365	2
366	5
367	6
368	6
369	7
370	4
371	3
372	3
373	10
374	8
375	7
376	5
377	4
378	5
379	3
380	6
381	3
382	4
383	6
384	4
385	1
386	5
387	3
388	5
389	8
390	3
391	5
392	3
393	5
394	2
395	7
396	3
397	6
398	7
399	3
400	3
401	3
402	5
403	3
404	7
405	6
406	1
407	7
408	7
409	6
410	5
411	7
412	2
413	6
414	9
415	3
416	3
417	4
418	5
419	3
420	6
421	1
422	4
423	5
424	7
425	4
426	5
427	7
428	4
429	7
430	6
431	11
432	7
433	3
434	5
435	7
436	4
437	3
438	6
439	4
440	2
441	5
442	6
443	5
444	3
445	3
446	3
447	5
448	4
449	7
450	2
451	5
452	3
453	6
454	2
455	4
456	5
457	8
458	4
459	5
460	3
461	3
462	1
463	4
464	3
465	4
466	1
467	2
468	6
469	4
470	4
471	3
472	3
473	4
474	2
475	2
476	4
477	1
478	6
479	4
480	1
481	5
482	1
483	4
484	2
486	4
487	2
488	5
489	4
490	6
491	3
492	5
493	8
494	2
495	3
496	1
497	2
498	5
499	1
500	4
501	6
502	4
503	1
504	1
505	4
506	4
507	1
508	4
509	2
510	3
511	2
513	6
514	2
515	2
516	6
517	7
518	2
519	3
520	1
521	6
522	2
523	5
524	1
525	2
526	2
527	2
528	1
529	3
530	2
531	2
532	1
533	4
534	8
535	1
536	2
537	4
538	2
539	2
540	1
541	1
542	1
543	1
544	2
545	3
546	1
547	2
548	4
549	2
550	1
552	2
553	3
554	3
555	3
556	5
557	4
558	2
559	4
561	3
562	4
563	1
564	4
567	3
569	2
570	2
571	1
572	3
573	1
574	6
575	3
576	2
577	2
578	2
579	1
580	4
581	4
582	3
583	1
584	3
585	3
586	1
587	2
588	2
589	2
590	1
591	3
592	1
593	2
594	6
595	4
596	6
597	3
598	2
599	2
601	5
603	2
604	5
605	2
606	2
607	3
608	2
609	5
610	1
611	2
612	2
613	4
616	1
617	2
618	1
619	1
620	2
621	2
622	1
623	1
624	3
625	5
626	1
627	2
628	4
630	1
631	2
632	2
633	4
635	4
636	3
637	2
638	2
639	4
640	3
641	1
642	2
644	1
645	4
646	4
647	2
648	1
649	2
650	2
651	2
652	3
653	2
654	4
655	1
656	2
657	1
658	2
659	3
660	1
661	5
662	1
663	3
664	3
665	1
666	2
667	1
668	1
671	1
672	2
673	1
674	1
675	1
677	2
678	2
679	3
680	1
682	4
684	2
685	1
686	1
687	1
688	2
689	1
690	1
691	1
694	4
695	2
696	2
697	1
699	2
700	2
702	1
703	1
704	1
705	1
706	1
707	1
712	1
713	3
714	1
715	1
716	2
717	1
718	1
719	2
720	2
721	3
722	2
723	1
724	2
725	1
729	1
730	1
731	1
732	3
733	2
735	1
736	2
737	1
738	1
739	1
741	1
744	3
746	2
747	1
748	1
749	1
750	1
751	1
752	1
753	2
754	1
755	1
756	1
760	1
761	1
763	3
764	1
765	1
766	1
767	1
768	2
770	3
771	4
772	1
774	1
775	1
777	1
778	1
779	1
780	2
782	1
784	2
789	1
790	1
792	2
793	4
794	1
797	1
800	1
802	2
803	1
804	1
806	2
807	1
808	3
809	2
810	3
811	3
813	1
814	1
815	1
816	1
817	1
818	2
819	1
820	1
821	1
822	1
823	1
824	3
828	1
829	2
832	1
835	3
836	1
837	2
838	1
840	1
844	1
847	2
849	1
850	1
852	1
854	3
858	1
859	2
861	1
862	1
863	3
864	1
866	1
867	2
869	2
871	1
872	1
873	2
876	4
878	2
880	1
882	2
885	2
886	1
887	1
889	1
890	1
891	1
892	1
895	2
897	1
898	2
899	1
902	2
903	1
905	2
907	2
909	2
910	1
911	1
912	1
917	1
918	1
921	1
922	2
923	2
927	2
928	1
930	1
932	1
937	1
939	1
940	1
944	1
945	1
947	1
951	2
952	1
953	2
956	2
959	1
961	1
963	1
964	2
965	2
967	1
970	1
971	1
981	2
982	2
984	2
985	1
986	1
989	1
990	1
992	1
993	2
994	1
995	1
996	1
997	1
998	1
1000	2
1002	1
1003	1
1004	2
1005	1
1012	1
1013	1
1017	1
1021	2
1022	1
1025	1
1027	1
1030	1
1031	2
1033	2
1034	1
1036	1
1037	1
1038	2
1040	1
1041	1
1046	1
1051	1
1053	1
1056	1
1058	1
1061	1
1062	1
1063	2
1065	1
1067	1
1069	2
1070	1
1072	1
1073	1
1074	1
1075	1
1077	1
1079	1
1082	1
1083	1
1084	2
1086	1
1093	1
1095	2
1096	1
1099	1
1101	1
1109	1
1110	1
1115	1
1122	1
1124	2
1125	1
1127	1
1130	1
1132	1
1139	2
1140	2
1142	1
1143	1
1149	1
1151	1
1153	2
1156	1
1158	1
1162	1
1163	1
1164	1
1167	2
1168	1
1170	1
1177	1
1180	1
1183	1
1185	1
1189	2
1200	1
1202	1
1203	1
1204	1
1206	1
1212	2
1215	1
1216	1
1217	2
1218	1
1219	1
1224	1
1225	1
1228	1
1231	1
1238	1
1243	1
1246	1
1249	1
1259	1
1264	1
1269	2
1274	1
1283	1
1284	1
1286	1
1288	1
1289	1
1295	1
1305	1
1313	1
1317	1
1319	1
1326	1
1328	1
1341	1
1345	1
1350	1
1351	1
1354	2
1359	1
1365	1
1371	1
1375	1
1386	1
1388	1
1394	1
1422	1
1426	1
1427	1
1433	1
1435	1
1438	1
1439	1
1445	1
1452	1
1456	2
1462	1
1465	1
1472	1
1476	1
1478	1
1485	1
1486	1
1487	1
1492	1
1493	1
1499	1
1516	1
1533	1
1535	1
1536	1
1549	1
1550	1
1556	1
1558	1
1569	1
1588	2
1609	1
1629	1
1630	1
1644	2
1649	1
1652	1
1660	1
1686	2
1689	1
1692	1
1703	1
1712	1
1713	1
1719	1
1727	1
1730	1
1738	1
1741	1
1757	1
1759	1
1762	1
1772	1
1792	1
1824	1
1831	1
1833	1
1847	1
1855	1
1870	1
1878	1
1880	1
1900	1
1920	1
1924	1
1927	1
1934	2
1941	1
1970	1
1977	1
1987	1
1991	1
1999	1
2013	1
2016	1
2034	1
2036	1
2046	1
2047	1
2050	1
2054	1
2072	1
2076	1
2103	1
2108	1
2152	1
2154	1
2173	1
2176	1
2184	1
2186	1
2188	1
2236	1
2273	1
2282	1
2296	1
2414	1
2417	1
2443	1
2489	1
2548	1
2558	1
2570	2
2586	1
2590	1
2676	1
2690	1
2697	2
2714	1
2732	1
2736	1
2859	1
2900	1
2906	1
2935	1
2941	1
2965	1
2969	1
3049	1
3104	1
3109	1
3177	1
3432	1
3478	1
3536	1
3644	1
3663	1
4052	1
4500	1
4782	1
5341	1
5798	1
5871	1
5948	1
6559	1
6926	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
.............................................._........_..............................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:02:58. Running peak memory: 3.993GB.  
  PID: 340488;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort_dups.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam) > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_preseq_counts.txt` (340652)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 3.993GB.  
  PID: 340652;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_preseq_plot` (340695)
<pre>
Processing Jurkat_ChRO-seq_2
INFO: Found real counts for Jurkat_ChRO-seq_2 - Total (M): 33.989659 Unique (M): 30.4684

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.993GB.  
  PID: 340695;	Command: Rscript;	Return code: 0;	Memory used: 0.205GB

> `Library complexity`	QC_hg38/Jurkat_ChRO-seq_2_preseq_plot.pdf	Library complexity	QC_hg38/Jurkat_ChRO-seq_2_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.6458	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-18 12:44:23) elapsed: 392.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam -c 8 -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_bamQC.tsv` (340714)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/tmp_Jurkat_ChRO-seq_2_sort_bi4oecyh'
Processing with 8 cores...
Discarding 98 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr17_KI270730v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270509v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270515v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270750v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1']
Keeping 97 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270510v1', 'chrUn_KI270518v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270366v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 3.993GB.  
  PID: 340714;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.085GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_bamQC.tsv`

> `NRF`	0.53	PEPPRO	_RES_

> `PBC1`	0.8	PEPPRO	_RES_

> `PBC2`	8.39	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_unmap.bam`  

> `samtools view -b -@ 8 -f 4  /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_unmap.bam` (340990)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.993GB.  
  PID: 340990;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_temp.bam`

> `Unmapped_reads`	2435871	PEPPRO	_RES_

### Split BAM by strand (02-18 12:45:24) elapsed: 61.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_plus.bam` (341028)
<pre>
</pre>
Command completed. Elapsed time: 0:01:43. Running peak memory: 3.993GB.  
  PID: 341028;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_minus.bam` (341123)
<pre>
</pre>
Command completed. Elapsed time: 0:01:42. Running peak memory: 3.993GB.  
  PID: 341123;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-18 12:48:49) elapsed: 205.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (341212)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.993GB.  
  PID: 341212;	Command: sed;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_plus_TssEnrichment.txt` (341213)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.993GB.  
  PID: 341213;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.438GB


> `TSS_Plus_Score`	68.4	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_minus_TssEnrichment.txt` (341256)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 3.993GB.  
  PID: 341256;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.51GB


> `TSS_Minus_Score`	23.2	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_minus_TssEnrichment.txt` (341293)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.993GB.  
  PID: 341293;	Command: Rscript;	Return code: 0;	Memory used: 0.255GB

> `TSS enrichment`	QC_hg38/Jurkat_ChRO-seq_2_TSSenrichment.pdf	TSS enrichment	QC_hg38/Jurkat_ChRO-seq_2_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt` (341312,341313,341314,341315)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.993GB.  
  PID: 341312;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 341314;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 341313;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 341315;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_keep.txt` (341317)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.993GB.  
  PID: 341317;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-18 12:49:43) elapsed: 54.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_ensembl_tss.bed` (341319,341320)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.993GB.  
  PID: 341319;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 341320;	Command: bedtools;	Return code: 0;	Memory used: 0.097GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed` (341324,341325)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.993GB.  
  PID: 341324;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 341325;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_TSS_density.bed` (341327,341328,341329,341330)
<pre>
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 3.993GB.  
  PID: 341330;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 341327;	Command: bedtools;	Return code: 0;	Memory used: 0.029GB  
  PID: 341329;	Command: sort;	Return code: 0;	Memory used: 0.009GB  
  PID: 341328;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_gene_body_density.bed` (341611,341612,341613)
<pre>
</pre>
Command completed. Elapsed time: 0:01:03. Running peak memory: 3.993GB.  
  PID: 341612;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 341611;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB  
  PID: 341613;	Command: sort;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_pause_index.bed` (341669,341670,341671)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.993GB.  
  PID: 341669;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 341671;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 341670;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	78.33	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_pause_index.bed` (341677)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.993GB.  
  PID: 341677;	Command: Rscript;	Return code: 0;	Memory used: 0.245GB

> `Pause index`	QC_hg38/Jurkat_ChRO-seq_2_pause_index.pdf	Pause index	QC_hg38/Jurkat_ChRO-seq_2_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_pause_index.bed` (341696)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.993GB.  
  PID: 341696;	Command: pigz;	Return code: 0;	Memory used: 0.001GB


### Calculate Fraction of Reads in pre-mature mRNA (02-18 12:51:43) elapsed: 120.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_plus.bam`
34004054 12298505

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_minus.bam`
34004054 11929067

> `Minus_FRiP`	0.35	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_gene_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/signal_hg38/Jurkat_ChRO-seq_2_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_gene_sort.bed` (341779,341780)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.993GB.  
  PID: 341779;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 341780;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/signal_hg38/Jurkat_ChRO-seq_2_gene_coverage.bed` (341782)
<pre>
</pre>
Command completed. Elapsed time: 0:01:05. Running peak memory: 3.993GB.  
  PID: 341782;	Command: bedtools;	Return code: 0;	Memory used: 0.101GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/raw/hg38_annotations.bed.gz` (341841)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.993GB.  
  PID: 341841;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/raw/hg38_annotations.bed` (341842)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.993GB.  
  PID: 341842;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate fraction and proportion of reads in features (FRiF/PRiF) (02-18 12:53:50) elapsed: 127.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/raw/hg38_annotations.bed` (341850)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.993GB.  
  PID: 341850;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/3_UTR"` (341853)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.993GB.  
  PID: 341853;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/3_UTR_sort.bed` (341854,341855,341857,341858)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.993GB.  
  PID: 341854;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 341855;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 341858;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB  
  PID: 341857;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_3_UTR_plus_coverage.bed` (341861)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.993GB.  
  PID: 341861;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_3_UTR_minus_coverage.bed` (341881)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.993GB.  
  PID: 341881;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/5_UTR"` (341899)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.993GB.  
  PID: 341899;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/5_UTR_sort.bed` (341900,341901,341902,341903)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.993GB.  
  PID: 341900;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 341901;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 341903;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 341902;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_5_UTR_plus_coverage.bed` (341906)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 3.993GB.  
  PID: 341906;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_5_UTR_minus_coverage.bed` (341927)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 3.993GB.  
  PID: 341927;	Command: bedtools;	Return code: 0;	Memory used: 0.029GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Enhancer_sort.bed` (342169,342170,342171,342172)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.993GB.  
  PID: 342169;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 342170;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 342172;	Command: bedtools;	Return code: 0;	Memory used: 0.048GB  
  PID: 342171;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Enhancer_plus_coverage.bed` (342174)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.993GB.  
  PID: 342174;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Enhancer_minus_coverage.bed` (342192)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.993GB.  
  PID: 342192;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Exon_sort.bed` (342209,342210,342211,342212)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.993GB.  
  PID: 342209;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 342210;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 342212;	Command: bedtools;	Return code: 0;	Memory used: 0.174GB  
  PID: 342211;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Exon_plus_coverage.bed` (342217)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.993GB.  
  PID: 342217;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Exon_minus_coverage.bed` (342244)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.993GB.  
  PID: 342244;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Intron_sort.bed` (342267,342268,342269,342270)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.993GB.  
  PID: 342267;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 342270;	Command: bedtools;	Return code: 0;	Memory used: 0.084GB  
  PID: 342268;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 342269;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Intron_plus_coverage.bed` (342274)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.993GB.  
  PID: 342274;	Command: bedtools;	Return code: 0;	Memory used: 0.067GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Intron_minus_coverage.bed` (342300)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 3.993GB.  
  PID: 342300;	Command: bedtools;	Return code: 0;	Memory used: 0.072GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter_sort.bed` (342323,342324,342325,342326)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.993GB.  
  PID: 342323;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 342324;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 342326;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 342325;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Promoter_plus_coverage.bed` (342328)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.993GB.  
  PID: 342328;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Promoter_minus_coverage.bed` (342354)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.993GB.  
  PID: 342354;	Command: bedtools;	Return code: 0;	Memory used: 0.047GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter_Flanking_Region"` (342377)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.993GB.  
  PID: 342377;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed` (342378,342379,342380,342381)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.993GB.  
  PID: 342378;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 342380;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 342379;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 342381;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed` (342385)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.993GB.  
  PID: 342385;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Promoter_Flanking_Region_minus_coverage.bed` (342403)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.993GB.  
  PID: 342403;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB


### Plot FRiF/PRiF (02-18 12:59:22) elapsed: 332.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s Jurkat_ChRO-seq_2 -z 3099922541 -n 17068847 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_frif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed` (342436)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.993GB.  
  PID: 342436;	Command: Rscript;	Return code: 0;	Memory used: 0.496GB

> `FRiF`	QC_hg38/Jurkat_ChRO-seq_2_frif.pdf	FRiF	QC_hg38/Jurkat_ChRO-seq_2_frif.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_prif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s Jurkat_ChRO-seq_2 -z 3099922541 -n 17068847 -y prif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_prif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed` (342477)
<pre>
Cumulative prif plot completed!

</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.993GB.  
  PID: 342477;	Command: Rscript;	Return code: 0;	Memory used: 0.502GB

> `PRiF`	QC_hg38/Jurkat_ChRO-seq_2_prif.pdf	PRiF	QC_hg38/Jurkat_ChRO-seq_2_prif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-18 13:00:26) elapsed: 64.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_exons_sort.bed` (342759,342760)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.993GB.  
  PID: 342760;	Command: bedtools;	Return code: 0;	Memory used: 0.084GB  
  PID: 342759;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_introns_sort.bed` (342766,342767,342768)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.993GB.  
  PID: 342766;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 342768;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 342767;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_exons_coverage.bed` (342775)
<pre>
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 3.993GB.  
  PID: 342775;	Command: bedtools;	Return code: 0;	Memory used: 0.032GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_introns_coverage.bed` (342826)
<pre>
</pre>
Command completed. Elapsed time: 0:00:51. Running peak memory: 3.993GB.  
  PID: 342826;	Command: bedtools;	Return code: 0;	Memory used: 0.083GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/34.004054)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_exons_rpkm.bed` (342872,342873,342874)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.993GB.  
  PID: 342872;	Command: awk;	Return code: 0;	Memory used: 0.009GB  
  PID: 342874;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 342873;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/34.004054)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_introns_rpkm.bed` (342877,342878,342879)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.993GB.  
  PID: 342877;	Command: awk;	Return code: 0;	Memory used: 0.01GB  
  PID: 342879;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 342878;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_exon_intron_ratios.bed` (342881,342882,342883)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.993GB.  
  PID: 342881;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 342883;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 342882;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.2	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_exon_intron_ratios.bed --annotate` (342890)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.993GB.  
  PID: 342890;	Command: Rscript;	Return code: 0;	Memory used: 0.238GB

> `mRNA contamination`	QC_hg38/Jurkat_ChRO-seq_2_mRNA_contamination.pdf	mRNA contamination	QC_hg38/Jurkat_ChRO-seq_2_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/QC_hg38/Jurkat_ChRO-seq_2_exon_intron_ratios.bed` (342909)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.993GB.  
  PID: 342909;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (02-18 13:02:18) elapsed: 112.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/signal_hg38/Jurkat_ChRO-seq_2_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/signal_hg38/Jurkat_ChRO-seq_2_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_plus.bam` (342917)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.993GB.  
  PID: 342917;	Command: samtools;	Return code: 0;	Memory used: 0.007GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/signal_hg38/Jurkat_ChRO-seq_2_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/signal_hg38/Jurkat_ChRO-seq_2_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (342934)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_plus.bam'
Temporary files will be stored in: 'tmp_Jurkat_ChRO-seq_2_plus_cuttrace_p79dowbp'
Processing with 2 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 114 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270509v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2']
Keeping 81 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270510v1', 'chrUn_KI270518v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270507v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270591v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 81 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/signal_hg38/Jurkat_ChRO-seq_2_plus_exact_body_0-mer.bw'
Merging 81 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/signal_hg38/Jurkat_ChRO-seq_2_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:10:37. Running peak memory: 3.993GB.  
  PID: 342934;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.325GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/signal_hg38/Jurkat_ChRO-seq_2_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/signal_hg38/Jurkat_ChRO-seq_2_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_minus.bam` (344740)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.993GB.  
  PID: 344740;	Command: samtools;	Return code: 0;	Memory used: 0.007GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/signal_hg38/Jurkat_ChRO-seq_2_minus_exact_body_0-mer.bw -p 5 --variable-step --tail-edge` (344753)
<pre>
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/aligned_hg38/Jurkat_ChRO-seq_2_minus.bam'
Temporary files will be stored in: 'tmp_Jurkat_ChRO-seq_2_minus_cuttrace_hk7uv8oy'
Processing with 5 cores...
stdin is empty of data
Discarding 111 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270719v1_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr17_KI270730v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270515v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1']
Keeping 84 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270510v1', 'chrUn_KI270519v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270538v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270330v1', 'chrUn_KI270366v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/Jurkat_ChRO-seq_2/signal_hg38/Jurkat_ChRO-seq_2_minus_exact_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 3.993GB.  
  PID: 344753;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.306GB

Starting cleanup: 63 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  2:04:52
*  Total elapsed time (all runs):  3:10:43
*         Peak memory (this run):  3.9933 GB
*        Pipeline completed time: 2020-02-18 13:13:50
