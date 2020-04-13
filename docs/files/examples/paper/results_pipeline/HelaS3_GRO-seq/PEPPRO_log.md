### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name HelaS3_GRO-seq --genome hg38 --input /project/shefflab/data/sra_fastq/SRR1693611.fastq.gz /project/shefflab/data/sra_fastq/SRR1693612.fastq.gz --single-or-paired single --protocol GRO --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/dev4/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-ba26-18
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/
*  Pipeline started at:   (02-27 09:24:33) elapsed: 1.0 _TIME_

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
*              `cores`:  `8`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `False`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data/sra_fastq/SRR1693611.fastq.gz', '/project/shefflab/data/sra_fastq/SRR1693612.fastq.gz']`
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
*           `protocol`:  `GRO`
*            `recover`:  `False`
*        `sample_name`:  `HelaS3_GRO-seq`
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

Local input file: /project/shefflab/data/sra_fastq/SRR1693611.fastq.gz

> `File_mb`	3550.33	PEPPRO	_RES_

> `Read_type`	single	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected GRO input

### Merge/link and fastq conversion:  (02-27 09:24:34) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz`  

> `cat /project/shefflab/data/sra_fastq/SRR1693611.fastq.gz /project/shefflab/data/sra_fastq/SRR1693612.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz` (13720)
<pre>
</pre>
Command completed. Elapsed time: 0:00:56. Running peak memory: 0.003GB.  
  PID: 13720;	Command: cat;	Return code: 0;	Memory used: 0.003GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1.fastq` (13970)
<pre>
</pre>
Command completed. Elapsed time: 0:01:58. Running peak memory: 0.003GB.  
  PID: 13970;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	72718941	PEPPRO	_RES_

> `Fastq_reads`	72718941	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz']

### FASTQ processing:  (02-27 09:28:43) elapsed: 250.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_processed.fastq`  

> `(cutadapt -j 8 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt` (14173)
<pre>
</pre>
Command completed. Elapsed time: 0:02:22. Running peak memory: 0.292GB.  
  PID: 14173;	Command: cutadapt;	Return code: 0;	Memory used: 0.292GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_noadap.fastq | seqtk seq -L 2 - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_processed.fastq` (14548,14549)
<pre>
</pre>
Command completed. Elapsed time: 0:01:21. Running peak memory: 0.292GB.  
  PID: 14549;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB  
  PID: 14548;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	22179715.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Uninformative_adapter_reads`	0.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	0.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/fastqc/HelaS3_GRO-seq_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (14714)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.292GB.  
  PID: 14714;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	72718941	PEPPRO	_RES_

> `Trim_loss_rate`	0.0	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_processed.fastq` (14735)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of HelaS3_GRO-seq_R1_processed.fastq
Approx 5% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 10% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 15% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 20% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 25% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 30% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 35% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 40% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 45% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 50% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 55% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 60% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 65% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 70% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 75% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 80% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 85% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 90% complete for HelaS3_GRO-seq_R1_processed.fastq
Approx 95% complete for HelaS3_GRO-seq_R1_processed.fastq
Analysis complete for HelaS3_GRO-seq_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:04:33. Running peak memory: 0.292GB.  
  PID: 14735;	Command: fastqc;	Return code: 0;	Memory used: 0.207GB

> `FastQC report r1`	fastqc/HelaS3_GRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/fastq/processed_R1.flag` (15189)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.292GB.  
  PID: 15189;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (02-27 09:38:02) elapsed: 558.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt` (15190)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.292GB.  
  PID: 15190;	Command: Rscript;	Return code: 0;	Memory used: 0.214GB

> `Adapter insertion distribution`	cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	0	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-27 09:38:07) elapsed: 6.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk 'END {print $1}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 10 && $1 >= 14){degradedSum += $2}; ($1 >= 4 && $1 <= 0){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.0	PEPPRO	_RES_

### Prealignments (02-27 09:38:08) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-27 09:38:08) elapsed: 0.0 _TIME_


> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id HelaS3_GRO-seq -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
72718941 reads; of these:
  72718941 (100.00%) were unpaired; of these:
    31839185 (43.78%) aligned 0 times
    40879756 (56.22%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
56.22% overall alignment rate

> `Aligned_reads_human_rDNA`	40879756.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	56.22	PEPPRO	_RES_

### Map to genome (02-27 09:45:02) elapsed: 415.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id HelaS3_GRO-seq -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/tmpfpcddf_3 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam` (16026,16036,16038)
<pre>
31839185 reads; of these:
  31839185 (100.00%) were unpaired; of these:
    10809494 (33.95%) aligned 0 times
    13023651 (40.90%) aligned exactly 1 time
    8006040 (25.15%) aligned >1 times
66.05% overall alignment rate
[bam_sort_core] merging from 8 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:19:19. Running peak memory: 3.593GB.  
  PID: 16026;	Command: bowtie2;	Return code: 0;	Memory used: 3.593GB  
  PID: 16036;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 16038;	Command: samtools;	Return code: 0;	Memory used: 0.889GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_fail_qc.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam` (17762)
<pre>
</pre>
Command completed. Elapsed time: 0:01:45. Running peak memory: 3.593GB.  
  PID: 17762;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	21029691	PEPPRO	_RES_

> `QC_filtered_reads`	6163527	PEPPRO	_RES_

> `Aligned_reads`	14866164	PEPPRO	_RES_

> `Alignment_rate`	20.44	PEPPRO	_RES_

> `Total_efficiency`	20.44	PEPPRO	_RES_

> `Read_depth`	2.34	PEPPRO	_RES_

### Compress all unmapped read files (02-27 10:17:27) elapsed: 1945.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq` (19084)
<pre>
</pre>
Command completed. Elapsed time: 0:01:00. Running peak memory: 3.593GB.  
  PID: 19084;	Command: pigz;	Return code: 0;	Memory used: 0.007GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam` (19152)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 3.593GB.  
  PID: 19152;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	1906316	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam` (19192)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.593GB.  
  PID: 19192;	Command: samtools;	Return code: 0;	Memory used: 0.007GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_noMT.bam` (19204,19205,19206,19207)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.593GB.  
  PID: 19205;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 19204;	Command: samtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 19206;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 19207;	Command: xargs;	Return code: 0;	Memory used: 0.048GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_noMT.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam` (19234)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 19234;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam` (19235)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.593GB.  
  PID: 19235;	Command: samtools;	Return code: 0;	Memory used: 0.007GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	51	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-27 10:20:11) elapsed: 164.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -c 8 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv` (19500)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/tmp_HelaS3_GRO-seq_sort_uuz40wm0'
Processing with 8 cores...
Discarding 117 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr22_KI270735v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 78 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270538v1', 'chrUn_KI270582v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000224v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 3.593GB.  
  PID: 19500;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.744GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv`

> `NRF`	0.72	PEPPRO	_RES_

> `PBC1`	0.91	PEPPRO	_RES_

> `PBC2`	13.55	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_unmap.bam`  

> `samtools view -b -@ 8 -f 4  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_unmap.bam` (19542)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.593GB.  
  PID: 19542;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam`

> `Unmapped_reads`	10809494	PEPPRO	_RES_

### Split BAM by strand (02-27 10:20:51) elapsed: 40.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam` (19582)
<pre>
</pre>
Command completed. Elapsed time: 0:00:58. Running peak memory: 3.593GB.  
  PID: 19582;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam` (19635)
<pre>
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 3.593GB.  
  PID: 19635;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-27 10:22:43) elapsed: 112.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (19694)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 19694;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_plus_TssEnrichment.txt` (19695)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.593GB.  
  PID: 19695;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.411GB


> `TSS_coding_score`	7.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_minus_TssEnrichment.txt` (19722)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.593GB.  
  PID: 19722;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.558GB


> `TSS_non-coding_score`	4.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_minus_TssEnrichment.txt` (19746)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.593GB.  
  PID: 19746;	Command: Rscript;	Return code: 0;	Memory used: 0.215GB

> `TSS enrichment`	QC_hg38/HelaS3_GRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/HelaS3_GRO-seq_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt` (19766,19767,19768,19769)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 19766;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 19768;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 19767;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 19769;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt` (19771)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 19771;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-27 10:23:05) elapsed: 22.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_tss.bed` (19773,19774)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.593GB.  
  PID: 19773;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 19774;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed` (19778,19779)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 19778;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 19779;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSS_density.bed` (19781,19782,19783,19784)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 3.593GB.  
  PID: 19781;	Command: bedtools;	Return code: 0;	Memory used: 0.032GB  
  PID: 19783;	Command: sort;	Return code: 0;	Memory used: 0.01GB  
  PID: 19782;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 19784;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_gene_body_density.bed` (19804,19805,19806)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 3.593GB.  
  PID: 19805;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 19804;	Command: bedtools;	Return code: 0;	Memory used: 0.137GB  
  PID: 19806;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed` (19835,19836,19837)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 19835;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 19837;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 19836;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	6.18	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed` (19842)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.593GB.  
  PID: 19842;	Command: Rscript;	Return code: 0;	Memory used: 0.215GB

> `Pause index`	QC_hg38/HelaS3_GRO-seq_pause_index.pdf	Pause index	QC_hg38/HelaS3_GRO-seq_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed` (19862)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 19862;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (02-27 10:24:05) elapsed: 60.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam`
14866164 5457524

> `Plus_FRiP`	0.37	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam`
14866164 4922927

> `Minus_FRiP`	0.33	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_gene_sort.bed` (19902,19903)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 19903;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 19902;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_gene_coverage.bed` (19906)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.593GB.  
  PID: 19906;	Command: bedtools;	Return code: 0;	Memory used: 0.104GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed.gz` (19933)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 19933;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed` (19934)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 19934;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (02-27 10:24:59) elapsed: 53.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed` (19943)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.593GB.  
  PID: 19943;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR"` (20021)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 20021;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR_sort.bed` (20047,20048,20053,20054)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 20047;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 20048;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 20054;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 20053;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_plus_coverage.bed` (20142)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.593GB.  
  PID: 20142;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_minus_coverage.bed` (20157)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.593GB.  
  PID: 20157;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR"` (20166)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 20166;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR_sort.bed` (20167,20168,20169,20170)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 20167;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 20168;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 20170;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 20169;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_plus_coverage.bed` (20173)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.593GB.  
  PID: 20173;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_minus_coverage.bed` (20183)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.593GB.  
  PID: 20183;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer_sort.bed` (20195,20196,20197,20198)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 20195;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 20196;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 20198;	Command: bedtools;	Return code: 0;	Memory used: 0.049GB  
  PID: 20197;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_plus_coverage.bed` (20202)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.593GB.  
  PID: 20202;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_minus_coverage.bed` (20212)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.593GB.  
  PID: 20212;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon_sort.bed` (20222,20223,20224,20225)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.593GB.  
  PID: 20222;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 20223;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 20225;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 20224;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_plus_coverage.bed` (20230)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.593GB.  
  PID: 20230;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_minus_coverage.bed` (20243)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.593GB.  
  PID: 20243;	Command: bedtools;	Return code: 0;	Memory used: 0.056GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron_sort.bed` (20255,20256,20257,20258)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.593GB.  
  PID: 20255;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 20257;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 20256;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 20258;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_plus_coverage.bed` (20261)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.593GB.  
  PID: 20261;	Command: bedtools;	Return code: 0;	Memory used: 0.04GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_minus_coverage.bed` (20277)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.593GB.  
  PID: 20277;	Command: bedtools;	Return code: 0;	Memory used: 0.049GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_sort.bed` (20288,20290,20291,20292)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 20288;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 20290;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 20292;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 20291;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_plus_coverage.bed` (20294)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.593GB.  
  PID: 20294;	Command: bedtools;	Return code: 0;	Memory used: 0.061GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_minus_coverage.bed` (20306)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.593GB.  
  PID: 20306;	Command: bedtools;	Return code: 0;	Memory used: 0.067GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region"` (20315)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 20315;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed` (20316,20317,20318,20319)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 20316;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 20318;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 20317;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 20319;	Command: bedtools;	Return code: 0;	Memory used: 0.05GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (20322)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.593GB.  
  PID: 20322;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_minus_coverage.bed` (20333)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.593GB.  
  PID: 20333;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB


### Plot cFRiF/FRiF (02-27 10:27:44) elapsed: 165.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s HelaS3_GRO-seq -z 3099922541 -n 7116450 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (20358)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.593GB.  
  PID: 20358;	Command: Rscript;	Return code: 0;	Memory used: 0.466GB

> `cFRiF`	QC_hg38/HelaS3_GRO-seq_cFRiF.pdf	cFRiF	QC_hg38/HelaS3_GRO-seq_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s HelaS3_GRO-seq -z 3099922541 -n 7116450 -y frif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (20397)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 3.593GB.  
  PID: 20397;	Command: Rscript;	Return code: 0;	Memory used: 0.466GB

> `FRiF`	QC_hg38/HelaS3_GRO-seq_FRiF.pdf	FRiF	QC_hg38/HelaS3_GRO-seq_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-27 10:28:46) elapsed: 62.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_exons_sort.bed` (20434,20435)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.593GB.  
  PID: 20435;	Command: bedtools;	Return code: 0;	Memory used: 0.088GB  
  PID: 20434;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_introns_sort.bed` (20443,20444,20445)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.593GB.  
  PID: 20444;	Command: bedtools;	Return code: 0;	Memory used: 0.078GB  
  PID: 20443;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 20445;	Command: bedtools;	Return code: 0;	Memory used: 0.092GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_coverage.bed` (20453)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 3.593GB.  
  PID: 20453;	Command: bedtools;	Return code: 0;	Memory used: 0.391GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_coverage.bed` (20474)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 3.593GB.  
  PID: 20474;	Command: bedtools;	Return code: 0;	Memory used: 0.05GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/14.866164)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_rpkm.bed` (20498,20499,20500)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 20498;	Command: awk;	Return code: 0;	Memory used: 0.006GB  
  PID: 20500;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 20499;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/14.866164)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_rpkm.bed` (20502,20503,20504)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 20502;	Command: awk;	Return code: 0;	Memory used: 0.005GB  
  PID: 20504;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 20503;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_rpkm.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed` (20507,20508,20509)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 20507;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 20509;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 20508;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	3.22	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed --annotate` (20515)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.593GB.  
  PID: 20515;	Command: Rscript;	Return code: 0;	Memory used: 0.235GB

> `mRNA contamination`	QC_hg38/HelaS3_GRO-seq_mRNA_contamination.pdf	mRNA contamination	QC_hg38/HelaS3_GRO-seq_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed` (20534)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 20534;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (02-27 10:29:54) elapsed: 67.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam` (20542)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.593GB.  
  PID: 20542;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_smooth_body_0-mer.bw -p 5 --variable-step` (20548)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam'
Temporary files will be stored in: 'tmp_HelaS3_GRO-seq_plus_cuttrace_24y_wcu1'
Processing with 2 cores...
Discarding 134 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr14_GL000225v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 61 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 61 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_exact_body_0-mer.bw'
Merging 61 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:42. Running peak memory: 3.593GB.  
  PID: 20548;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.567GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam` (22025)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.593GB.  
  PID: 22025;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_smooth_body_0-mer.bw -p 5 --variable-step` (22032)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam'
Temporary files will be stored in: 'tmp_HelaS3_GRO-seq_minus_cuttrace_sev7_mwr'
Processing with 2 cores...
stdin is empty of data
Discarding 135 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_GL000216v2', 'chrEBV']
Keeping 60 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270538v1', 'chrUn_KI270582v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270742v1', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 60 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_exact_body_0-mer.bw'
Merging 60 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:43. Running peak memory: 3.593GB.  
  PID: 22032;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.52GB

Starting cleanup: 52 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  1:23:01
*  Total elapsed time (all runs):  1:37:17
*         Peak memory (this run):  3.5925 GB
*        Pipeline completed time: 2020-02-27 10:47:33
