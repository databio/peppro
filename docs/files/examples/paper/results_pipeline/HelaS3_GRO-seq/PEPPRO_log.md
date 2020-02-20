### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name HelaS3_GRO-seq --genome hg38 --input /project/shefflab/data/sra_fastq/SRR1693611.fastq.gz /project/shefflab/data/sra_fastq/SRR1693612.fastq.gz --single-or-paired single --protocol GRO --max-len -1 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-aj38-13c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/
*  Pipeline started at:   (02-18 11:08:57) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data/sra_fastq/SRR1693611.fastq.gz', '/project/shefflab/data/sra_fastq/SRR1693612.fastq.gz']`
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

### Merge/link and fastq conversion:  (02-18 11:08:57) elapsed: 0.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz`  

> `cat /project/shefflab/data/sra_fastq/SRR1693611.fastq.gz /project/shefflab/data/sra_fastq/SRR1693612.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz` (116508)
<pre>
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 0.003GB.  
  PID: 116508;	Command: cat;	Return code: 0;	Memory used: 0.003GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1.fastq` (116603)
<pre>
</pre>
Command completed. Elapsed time: 0:04:01. Running peak memory: 0.003GB.  
  PID: 116603;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	72718941	PEPPRO	_RES_

> `Fastq_reads`	72718941	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz']

### FASTQ processing:  (02-18 11:15:30) elapsed: 393.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_processed.fastq`  

> `(cutadapt -j 8 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt` (117453)
<pre>
</pre>
Command completed. Elapsed time: 0:03:30. Running peak memory: 0.297GB.  
  PID: 117453;	Command: cutadapt;	Return code: 0;	Memory used: 0.297GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_noadap.fastq | seqtk seq -L 2 - > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_processed.fastq` (117710,117712)
<pre>
</pre>
Command completed. Elapsed time: 0:02:12. Running peak memory: 0.297GB.  
  PID: 117712;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB  
  PID: 117710;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	22179715.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_too_short`	0.0	PEPPRO	_RES_

> `Pct_reads_too_short`	0.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/fastqc/HelaS3_GRO-seq_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (118212)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.297GB.  
  PID: 118212;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	72718941	PEPPRO	_RES_

> `Trim_loss_rate`	0.0	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_processed.fastq` (118230)
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
Command completed. Elapsed time: 0:08:36. Running peak memory: 0.297GB.  
  PID: 118230;	Command: fastqc;	Return code: 0;	Memory used: 0.216GB

> `FastQC report r1`	fastqc/HelaS3_GRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_

### Plot adapter insertion distribution (02-18 11:31:37) elapsed: 966.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt` (119262)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.297GB.  
  PID: 119262;	Command: Rscript;	Return code: 0;	Memory used: 0.127GB

> `Adapter insertion distribution`	cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	0	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-18 11:31:43) elapsed: 7.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk 'END {print $1}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 10 && $1 >= 14){degradedSum += $2}; ($1 >= 4 && $1 <= 0){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.0	PEPPRO	_RES_

### Prealignments (02-18 11:31:44) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-18 11:31:44) elapsed: 0.0 _TIME_


> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id HelaS3_GRO-seq -U /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
72718941 reads; of these:
  72718941 (100.00%) were unpaired; of these:
    31839185 (43.78%) aligned 0 times
    40879756 (56.22%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
56.22% overall alignment rate

> `Aligned_reads_human_rDNA`	40879756.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	56.22	PEPPRO	_RES_

### Map to genome (02-18 11:38:58) elapsed: 434.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id HelaS3_GRO-seq -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/tmpd0x3clbl -o /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam` (120054,120055,120056)
<pre>
31839185 reads; of these:
  31839185 (100.00%) were unpaired; of these:
    10809494 (33.95%) aligned 0 times
    13023651 (40.90%) aligned exactly 1 time
    8006040 (25.15%) aligned >1 times
66.05% overall alignment rate
[bam_sort_core] merging from 8 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:20:04. Running peak memory: 3.593GB.  
  PID: 120054;	Command: bowtie2;	Return code: 0;	Memory used: 3.593GB  
  PID: 120055;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 120056;	Command: samtools;	Return code: 0;	Memory used: 0.874GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam` (122309)
<pre>
</pre>
Command completed. Elapsed time: 0:01:43. Running peak memory: 3.593GB.  
  PID: 122309;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	21029691	PEPPRO	_RES_

> `QC_filtered_reads`	6163527	PEPPRO	_RES_

> `Aligned_reads`	14866164	PEPPRO	_RES_

> `Alignment_rate`	20.44	PEPPRO	_RES_

> `Total_efficiency`	20.44	PEPPRO	_RES_

> `Read_depth`	2.34	PEPPRO	_RES_

### Compress all unmapped read files (02-18 12:10:09) elapsed: 1872.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq` (123777)
<pre>
</pre>
Command completed. Elapsed time: 0:01:04. Running peak memory: 3.593GB.  
  PID: 123777;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam` (123880)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 3.593GB.  
  PID: 123880;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	1906316	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam` (123911)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.593GB.  
  PID: 123911;	Command: samtools;	Return code: 0;	Memory used: 0.007GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_noMT.bam` (123928,123929,123930,123931)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 3.593GB.  
  PID: 123929;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 123928;	Command: samtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 123930;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 123931;	Command: xargs;	Return code: 0;	Memory used: 0.048GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam` (123958)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 123958;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam` (123959)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.593GB.  
  PID: 123959;	Command: samtools;	Return code: 0;	Memory used: 0.007GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	51	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-18 12:12:52) elapsed: 163.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -c 8 -o /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv` (124026)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/tmp_HelaS3_GRO-seq_sort_gsd8f38g'
Processing with 8 cores...
Discarding 117 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr22_KI270735v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 78 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270538v1', 'chrUn_KI270582v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000224v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.593GB.  
  PID: 124026;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.751GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv`

> `NRF`	0.72	PEPPRO	_RES_

> `PBC1`	0.91	PEPPRO	_RES_

> `PBC2`	13.55	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_unmap.bam`  

> `samtools view -b -@ 8 -f 4  /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_unmap.bam` (124065)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.593GB.  
  PID: 124065;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam`

> `Unmapped_reads`	10809494	PEPPRO	_RES_

### Split BAM by strand (02-18 12:13:31) elapsed: 38.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam` (124106)
<pre>
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 3.593GB.  
  PID: 124106;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam` (124171)
<pre>
</pre>
Command completed. Elapsed time: 0:00:51. Running peak memory: 3.593GB.  
  PID: 124171;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-18 12:15:15) elapsed: 105.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (124442)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 124442;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_plus_TssEnrichment.txt` (124443)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.593GB.  
  PID: 124443;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.471GB


> `TSS_Plus_Score`	7.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_minus_TssEnrichment.txt` (124469)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.593GB.  
  PID: 124469;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.595GB


> `TSS_Minus_Score`	4.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_minus_TssEnrichment.txt` (124495)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.593GB.  
  PID: 124495;	Command: Rscript;	Return code: 0;	Memory used: 0.125GB

> `TSS enrichment`	QC_hg38/HelaS3_GRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/HelaS3_GRO-seq_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt` (124517,124518,124519,124520)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 124517;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 124519;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 124518;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 124520;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt` (124522)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 124522;	Command: cut;	Return code: 0;	Memory used: 0.003GB


### Calculate Pause Index (PI) (02-18 12:15:35) elapsed: 20.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_tss.bed` (124524,124525)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.593GB.  
  PID: 124524;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 124525;	Command: bedtools;	Return code: 0;	Memory used: 0.097GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed` (124528,124529)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 124528;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 124529;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSS_density.bed` (124532,124533,124534,124535)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 3.593GB.  
  PID: 124532;	Command: bedtools;	Return code: 0;	Memory used: 0.056GB  
  PID: 124534;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 124533;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 124535;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_gene_body_density.bed` (124563,124564,124565)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 3.593GB.  
  PID: 124564;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 124563;	Command: bedtools;	Return code: 0;	Memory used: 0.239GB  
  PID: 124565;	Command: sort;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed` (124596,124597,124598)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 124596;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 124598;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 124597;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	6.18	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed` (124604)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.593GB.  
  PID: 124604;	Command: Rscript;	Return code: 0;	Memory used: 0.235GB

> `Pause index`	QC_hg38/HelaS3_GRO-seq_pause_index.pdf	Pause index	QC_hg38/HelaS3_GRO-seq_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed` (124626)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 124626;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate Fraction of Reads in pre-mature mRNA (02-18 12:16:27) elapsed: 52.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam`
14866164 5457524

> `Plus_FRiP`	0.37	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam`
14866164 4922927

> `Minus_FRiP`	0.33	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_gene_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_gene_sort.bed` (124684,124685)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 124684;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 124685;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_gene_coverage.bed` (124687)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 3.593GB.  
  PID: 124687;	Command: bedtools;	Return code: 0;	Memory used: 0.104GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed.gz` (124714)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 124714;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed` (124715)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 124715;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate fraction and proportion of reads in features (FRiF/PRiF) (02-18 12:17:19) elapsed: 52.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed` (124724)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.593GB.  
  PID: 124724;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR"` (124729)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 124729;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR_sort.bed` (124730,124731,124732,124733)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 124730;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 124731;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 124733;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 124732;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_plus_coverage.bed` (124736)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.593GB.  
  PID: 124736;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_minus_coverage.bed` (124758)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.593GB.  
  PID: 124758;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR"` (124767)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 124767;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR_sort.bed` (124768,124769,124770,124771)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 124768;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 124769;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 124771;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 124770;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_plus_coverage.bed` (124774)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.593GB.  
  PID: 124774;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_minus_coverage.bed` (124790)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.593GB.  
  PID: 124790;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer_sort.bed` (124831,124832,124833,124834)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 124831;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 124832;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 124834;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 124833;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_plus_coverage.bed` (124837)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.593GB.  
  PID: 124837;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_minus_coverage.bed` (124846)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.593GB.  
  PID: 124846;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon_sort.bed` (124855,124856,124857,124858)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.593GB.  
  PID: 124855;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 124856;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 124858;	Command: bedtools;	Return code: 0;	Memory used: 0.175GB  
  PID: 124857;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_plus_coverage.bed` (124862)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.593GB.  
  PID: 124862;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_minus_coverage.bed` (124882)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.593GB.  
  PID: 124882;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron_sort.bed` (124892,124893,124894,124895)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 124892;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 124894;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 124893;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 124895;	Command: bedtools;	Return code: 0;	Memory used: 0.085GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_plus_coverage.bed` (124898)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.593GB.  
  PID: 124898;	Command: bedtools;	Return code: 0;	Memory used: 0.041GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_minus_coverage.bed` (124910)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.593GB.  
  PID: 124910;	Command: bedtools;	Return code: 0;	Memory used: 0.05GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_sort.bed` (124920,124921,124922,124923)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 124920;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 124921;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 124923;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 124922;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_plus_coverage.bed` (124926)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.593GB.  
  PID: 124926;	Command: bedtools;	Return code: 0;	Memory used: 0.062GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_minus_coverage.bed` (124942)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.593GB.  
  PID: 124942;	Command: bedtools;	Return code: 0;	Memory used: 0.068GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region"` (124954)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 124954;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed` (124955,124956,124957,124958)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 124955;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 124957;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 124956;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 124958;	Command: bedtools;	Return code: 0;	Memory used: 0.049GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (124961)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.593GB.  
  PID: 124961;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_minus_coverage.bed` (124971)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.593GB.  
  PID: 124971;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB


### Plot FRiF/PRiF (02-18 12:19:51) elapsed: 152.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s HelaS3_GRO-seq -z 3099922541 -n 7116450 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_frif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (124992)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 3.593GB.  
  PID: 124992;	Command: Rscript;	Return code: 0;	Memory used: 0.487GB

> `FRiF`	QC_hg38/HelaS3_GRO-seq_frif.pdf	FRiF	QC_hg38/HelaS3_GRO-seq_frif.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_prif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s HelaS3_GRO-seq -z 3099922541 -n 7116450 -y prif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_prif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (125253)
<pre>
Cumulative prif plot completed!

</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 3.593GB.  
  PID: 125253;	Command: Rscript;	Return code: 0;	Memory used: 0.487GB

> `PRiF`	QC_hg38/HelaS3_GRO-seq_prif.pdf	PRiF	QC_hg38/HelaS3_GRO-seq_prif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-18 12:20:47) elapsed: 56.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_exons_sort.bed` (125293,125294)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.593GB.  
  PID: 125294;	Command: bedtools;	Return code: 0;	Memory used: 0.094GB  
  PID: 125293;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_introns_sort.bed` (125300,125301,125302)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.593GB.  
  PID: 125300;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 125302;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB  
  PID: 125301;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_coverage.bed` (125308)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.593GB.  
  PID: 125308;	Command: bedtools;	Return code: 0;	Memory used: 0.396GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_coverage.bed` (125326)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 3.593GB.  
  PID: 125326;	Command: bedtools;	Return code: 0;	Memory used: 0.091GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/14.866164)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_rpkm.bed` (125353,125354,125355)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 125353;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 125355;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 125354;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/14.866164)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_rpkm.bed` (125357,125358,125359)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.593GB.  
  PID: 125357;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 125359;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 125358;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed` (125362,125363,125364)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 125362;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 125364;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 125363;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	3.22	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed --annotate` (125370)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.593GB.  
  PID: 125370;	Command: Rscript;	Return code: 0;	Memory used: 0.324GB

> `mRNA contamination`	QC_hg38/HelaS3_GRO-seq_mRNA_contamination.pdf	mRNA contamination	QC_hg38/HelaS3_GRO-seq_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed` (125389)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.593GB.  
  PID: 125389;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (02-18 12:21:44) elapsed: 57.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam` (125398)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.593GB.  
  PID: 125398;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_smooth_body_0-mer.bw -p 5 --variable-step` (125403)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam'
Temporary files will be stored in: 'tmp_HelaS3_GRO-seq_plus_cuttrace_v8mg7ruj'
Processing with 2 cores...
Discarding 134 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr14_GL000225v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 61 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 61 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_exact_body_0-mer.bw'
Merging 61 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:51. Running peak memory: 3.593GB.  
  PID: 125403;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.076GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam` (126693)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.593GB.  
  PID: 126693;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_exact_body_0-mer.bw -p 5 --variable-step` (126698)
<pre>
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam'
Temporary files will be stored in: 'tmp_HelaS3_GRO-seq_minus_cuttrace_p1268usx'
Processing with 5 cores...
stdin is empty of data
Discarding 135 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_GL000216v2', 'chrEBV']
Keeping 60 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270538v1', 'chrUn_KI270582v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270742v1', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 60 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_exact_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.593GB.  
  PID: 126698;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.29GB

Starting cleanup: 52 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  1:21:02
*  Total elapsed time (all runs):  1:38:02
*         Peak memory (this run):  3.5925 GB
*        Pipeline completed time: 2020-02-18 12:29:59
