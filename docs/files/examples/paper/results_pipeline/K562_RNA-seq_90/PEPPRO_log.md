### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name K562_RNA-seq_90 --genome hg38 --input /project/shefflab/data/guertin/fastq/K562_90pct_RNArc_r2.fastq.gz --single-or-paired single --protocol PRO --max-len -1 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/results_pipeline -P 4 -M 16000`
*         Compute host:  udc-aj38-18c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/
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
*              `cores`:  `4`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `False`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data/guertin/fastq/K562_90pct_RNArc_r2.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_90`
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

Local input file: /project/shefflab/data/guertin/fastq/K562_90pct_RNArc_r2.fastq.gz

> `File_mb`	779.46	PEPPRO	_RES_

> `Read_type`	single	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-18 11:08:58) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/raw/K562_RNA-seq_90.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/K562_90pct_RNArc_r2.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/raw/K562_RNA-seq_90.fastq.gz` (105976)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.002GB.  
  PID: 105976;	Command: ln;	Return code: 0;	Memory used: 0.002GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/raw/K562_RNA-seq_90.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/fastq/K562_RNA-seq_90_R1.fastq`  

> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/raw/K562_RNA-seq_90.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/fastq/K562_RNA-seq_90_R1.fastq` (105977)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 0.003GB.  
  PID: 105977;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	10000000	PEPPRO	_RES_

> `Fastq_reads`	10000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/raw/K562_RNA-seq_90.fastq.gz']

### FASTQ processing:  (02-18 11:09:40) elapsed: 42.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/fastq/K562_RNA-seq_90_R1_processed.fastq`  

> `(cutadapt -j 4 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/fastq/K562_RNA-seq_90_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/fastq/K562_RNA-seq_90_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/cutadapt/K562_RNA-seq_90_R1_cutadapt.txt` (106045)
<pre>
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 0.167GB.  
  PID: 106045;	Command: cutadapt;	Return code: 0;	Memory used: 0.167GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/fastq/K562_RNA-seq_90_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/fastq/K562_RNA-seq_90_R1_processed.fastq` (106333,106334)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 0.167GB.  
  PID: 106333;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 106334;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/cutadapt/K562_RNA-seq_90_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	4024168.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/cutadapt/K562_RNA-seq_90_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/cutadapt/K562_RNA-seq_90_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/cutadapt/K562_RNA-seq_90_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_too_short`	24881.0	PEPPRO	_RES_

> `Pct_reads_too_short`	0.2488	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/fastqc/K562_RNA-seq_90_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (106382)
### Calculate the number of trimmed reads
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.167GB.  
  PID: 106382;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	9975119	PEPPRO	_RES_

> `Trim_loss_rate`	0.25	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/fastq/K562_RNA-seq_90_R1_processed.fastq` (106388)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of K562_RNA-seq_90_R1_processed.fastq
Approx 5% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 10% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 15% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 20% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 25% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 30% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 35% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 40% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 45% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 50% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 55% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 60% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 65% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 70% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 75% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 80% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 85% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 90% complete for K562_RNA-seq_90_R1_processed.fastq
Approx 95% complete for K562_RNA-seq_90_R1_processed.fastq
Analysis complete for K562_RNA-seq_90_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:45. Running peak memory: 0.196GB.  
  PID: 106388;	Command: fastqc;	Return code: 0;	Memory used: 0.196GB

> `FastQC report r1`	fastqc/K562_RNA-seq_90_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_

### Plot adapter insertion distribution (02-18 11:11:33) elapsed: 113.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/cutadapt/K562_RNA-seq_90_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/cutadapt/K562_RNA-seq_90_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/cutadapt` (106474)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 0.343GB.  
  PID: 106474;	Command: Rscript;	Return code: 0;	Memory used: 0.343GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_90_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_90_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/cutadapt/K562_RNA-seq_90_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-18 11:11:48) elapsed: 14.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/cutadapt/K562_RNA-seq_90_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/cutadapt/K562_RNA-seq_90_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/cutadapt/K562_RNA-seq_90_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/cutadapt/K562_RNA-seq_90_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/cutadapt/K562_RNA-seq_90_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.23	PEPPRO	_RES_

### Prealignments (02-18 11:11:48) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-18 11:11:48) elapsed: 0.0 _TIME_


> `(bowtie2 -p 4 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_90 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/fastq/K562_RNA-seq_90_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/prealignments/K562_RNA-seq_90_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
9975119 reads; of these:
  9975119 (100.00%) were unpaired; of these:
    9661759 (96.86%) aligned 0 times
    313360 (3.14%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
3.14% overall alignment rate

> `Aligned_reads_human_rDNA`	313360.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	3.14	PEPPRO	_RES_

### Map to genome (02-18 11:13:22) elapsed: 94.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam`  

> `bowtie2 -p 4 --very-sensitive -X 2000 --rg-id K562_RNA-seq_90 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/prealignments/K562_RNA-seq_90_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/tmpqfz36a3m -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_temp.bam` (106673,106674,106675)
<pre>
9661759 reads; of these:
  9661759 (100.00%) were unpaired; of these:
    1355949 (14.03%) aligned 0 times
    5045625 (52.22%) aligned exactly 1 time
    3260185 (33.74%) aligned >1 times
85.97% overall alignment rate
[bam_sort_core] merging from 2 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:12:58. Running peak memory: 3.53GB.  
  PID: 106673;	Command: bowtie2;	Return code: 0;	Memory used: 3.53GB  
  PID: 106674;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 106675;	Command: samtools;	Return code: 0;	Memory used: 0.865GB


> `samtools view -q 10 -b -@ 4 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam` (108297)
<pre>
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 3.53GB.  
  PID: 108297;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	8305810	PEPPRO	_RES_

> `QC_filtered_reads`	1759025	PEPPRO	_RES_

> `Aligned_reads`	6546785	PEPPRO	_RES_

> `Alignment_rate`	65.63	PEPPRO	_RES_

> `Total_efficiency`	65.47	PEPPRO	_RES_

> `Read_depth`	3.75	PEPPRO	_RES_

### Compress all unmapped read files (02-18 11:31:33) elapsed: 1090.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/prealignments/K562_RNA-seq_90_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/prealignments/K562_RNA-seq_90_human_rDNA_unmap.fq` (108917)
<pre>
</pre>
Command completed. Elapsed time: 0:00:59. Running peak memory: 3.53GB.  
  PID: 108917;	Command: pigz;	Return code: 0;	Memory used: 0.004GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_temp.bam` (108985)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.53GB.  
  PID: 108985;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	560951	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam` (108998)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.53GB.  
  PID: 108998;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_noMT.bam` (109005,109006,109007,109008)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 3.53GB.  
  PID: 109006;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 109005;	Command: samtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 109007;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 109008;	Command: xargs;	Return code: 0;	Memory used: 0.026GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam` (109044)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.53GB.  
  PID: 109044;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam` (109045)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.53GB.  
  PID: 109045;	Command: samtools;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-18 11:33:37) elapsed: 125.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam -c 4 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_bamQC.tsv` (109078)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/tmp_K562_RNA-seq_90_sort_1cgy7vwn'
Processing with 4 cores...
Discarding 122 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrEBV']
Keeping 73 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.53GB.  
  PID: 109078;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.317GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_bamQC.tsv`

> `NRF`	0.76	PEPPRO	_RES_

> `PBC1`	0.91	PEPPRO	_RES_

> `PBC2`	18.23	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_unmap.bam`  

> `samtools view -b -@ 4 -f 4  /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_unmap.bam` (109108)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.53GB.  
  PID: 109108;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools view -c -f 4 -@ 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_temp.bam`

> `Unmapped_reads`	1355949	PEPPRO	_RES_

### Split BAM by strand (02-18 11:34:07) elapsed: 30.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_plus.bam` (109131)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 3.53GB.  
  PID: 109131;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_minus.bam` (109171)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 3.53GB.  
  PID: 109171;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-18 11:35:13) elapsed: 66.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (109417)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.53GB.  
  PID: 109417;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/plus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_plus_TssEnrichment.txt` (109418)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.53GB.  
  PID: 109418;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.196GB


> `TSS_Plus_Score`	18.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/minus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_minus_TssEnrichment.txt` (109437)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.53GB.  
  PID: 109437;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.2GB


> `TSS_Minus_Score`	3.2	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_minus_TssEnrichment.txt` (109466)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.53GB.  
  PID: 109466;	Command: Rscript;	Return code: 0;	Memory used: 0.205GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_90_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_90_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt` (109484,109485,109486,109487)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.53GB.  
  PID: 109484;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 109486;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 109485;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 109487;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_keep.txt` (109489)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.53GB.  
  PID: 109489;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (02-18 11:35:35) elapsed: 22.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_ensembl_tss.bed` (109491,109492)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.53GB.  
  PID: 109491;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 109492;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_ensembl_gene_body.bed` (109496,109497)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.53GB.  
  PID: 109496;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 109497;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_TSS_density.bed` (109500,109501,109502,109503)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.53GB.  
  PID: 109502;	Command: sort;	Return code: 0;	Memory used: 0.008GB  
  PID: 109500;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 109503;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 109501;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_gene_body_density.bed` (109514,109515,109516)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.53GB.  
  PID: 109514;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB  
  PID: 109516;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 109515;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_pause_index.bed` (109538,109539,109540)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.53GB.  
  PID: 109538;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 109540;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 109539;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	12.14	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_pause_index.bed` (109545)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.53GB.  
  PID: 109545;	Command: Rscript;	Return code: 0;	Memory used: 0.236GB

> `Pause index`	QC_hg38/K562_RNA-seq_90_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_90_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_pause_index.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_pause_index.bed` (109562)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.53GB.  
  PID: 109562;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (02-18 11:36:10) elapsed: 35.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_plus.bam`
6546785 2597009

> `Plus_FRiP`	0.4	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_minus.bam`
6546785 2612686

> `Minus_FRiP`	0.4	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_gene_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/signal_hg38/K562_RNA-seq_90_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_gene_sort.bed` (109589,109590)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.53GB.  
  PID: 109589;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 109590;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/signal_hg38/K562_RNA-seq_90_gene_coverage.bed` (109593)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.53GB.  
  PID: 109593;	Command: bedtools;	Return code: 0;	Memory used: 0.058GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/raw/hg38_annotations.bed.gz` (109610)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.53GB.  
  PID: 109610;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/raw/hg38_annotations.bed` (109611)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.53GB.  
  PID: 109611;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate fraction and proportion of reads in features (FRiF/PRiF) (02-18 11:36:37) elapsed: 28.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/raw/hg38_annotations.bed` (109620)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.53GB.  
  PID: 109620;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/3_UTR"` (109622)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.53GB.  
  PID: 109622;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/3_UTR_sort.bed` (109623,109624,109625,109626)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.53GB.  
  PID: 109623;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 109624;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 109626;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 109625;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_3_UTR_plus_coverage.bed` (109629)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.53GB.  
  PID: 109629;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_3_UTR_minus_coverage.bed` (109637)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.53GB.  
  PID: 109637;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/5_UTR"` (109643)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.53GB.  
  PID: 109643;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/5_UTR_sort.bed` (109644,109645,109646,109647)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.53GB.  
  PID: 109644;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 109645;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 109647;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 109646;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_5_UTR_plus_coverage.bed` (109651)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.53GB.  
  PID: 109651;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_5_UTR_minus_coverage.bed` (109658)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.53GB.  
  PID: 109658;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Enhancer_sort.bed` (109664,109665,109666,109667)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.53GB.  
  PID: 109664;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 109665;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 109667;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB  
  PID: 109666;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Enhancer_plus_coverage.bed` (109670)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.53GB.  
  PID: 109670;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Enhancer_minus_coverage.bed` (109679)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.53GB.  
  PID: 109679;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Exon_sort.bed` (109699,109700,109701,109702)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.53GB.  
  PID: 109699;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 109700;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 109702;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 109701;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Exon_plus_coverage.bed` (109707)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.53GB.  
  PID: 109707;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Exon_minus_coverage.bed` (109715)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.53GB.  
  PID: 109715;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Intron_sort.bed` (109728,109729,109730,109731)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.53GB.  
  PID: 109728;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 109730;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 109729;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 109731;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Intron_plus_coverage.bed` (109734)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.53GB.  
  PID: 109734;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Intron_minus_coverage.bed` (109743)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.53GB.  
  PID: 109743;	Command: bedtools;	Return code: 0;	Memory used: 0.025GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter_sort.bed` (109759,109760,109761,109762)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.53GB.  
  PID: 109759;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 109760;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 109762;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 109761;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Promoter_plus_coverage.bed` (109765)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.53GB.  
  PID: 109765;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Promoter_minus_coverage.bed` (109776)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.53GB.  
  PID: 109776;	Command: bedtools;	Return code: 0;	Memory used: 0.023GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter_Flanking_Region"` (109786)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.53GB.  
  PID: 109786;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter_Flanking_Region_sort.bed` (109787,109788,109789,109790)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.53GB.  
  PID: 109787;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 109789;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 109788;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 109790;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Promoter_Flanking_Region_plus_coverage.bed` (109793)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.53GB.  
  PID: 109793;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Promoter_Flanking_Region_minus_coverage.bed` (109800)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.53GB.  
  PID: 109800;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB


### Plot FRiF/PRiF (02-18 11:38:14) elapsed: 96.0 _TIME_


> `samtools view -@ 4 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_90 -z 3099922541 -n 3050675 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_frif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Promoter_Flanking_Region_plus_coverage.bed` (109818)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 3.53GB.  
  PID: 109818;	Command: Rscript;	Return code: 0;	Memory used: 0.459GB

> `FRiF`	QC_hg38/K562_RNA-seq_90_frif.pdf	FRiF	QC_hg38/K562_RNA-seq_90_frif.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_prif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_90 -z 3099922541 -n 3050675 -y prif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_prif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_Promoter_Flanking_Region_plus_coverage.bed` (109878)
<pre>
Cumulative prif plot completed!

</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.53GB.  
  PID: 109878;	Command: Rscript;	Return code: 0;	Memory used: 0.459GB

> `PRiF`	QC_hg38/K562_RNA-seq_90_prif.pdf	PRiF	QC_hg38/K562_RNA-seq_90_prif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-18 11:39:12) elapsed: 58.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_exons_sort.bed` (109926,109927)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.53GB.  
  PID: 109927;	Command: bedtools;	Return code: 0;	Memory used: 0.087GB  
  PID: 109926;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_introns_sort.bed` (109933,109934,109935)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.53GB.  
  PID: 109933;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 109935;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 109934;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_exons_coverage.bed` (109941)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.53GB.  
  PID: 109941;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_introns_coverage.bed` (109965)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.53GB.  
  PID: 109965;	Command: bedtools;	Return code: 0;	Memory used: 0.04GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/6.546785)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_exons_rpkm.bed` (109988,109989,109990)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.53GB.  
  PID: 109988;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 109990;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 109989;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/6.546785)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_introns_rpkm.bed` (109992,109993,109994)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.53GB.  
  PID: 109992;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 109994;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 109993;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_exon_intron_ratios.bed` (109997,109998,109999)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.53GB.  
  PID: 109997;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 109999;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 109998;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	10.06	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_exon_intron_ratios.bed --annotate` (110005)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.53GB.  
  PID: 110005;	Command: Rscript;	Return code: 0;	Memory used: 0.272GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_90_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_90_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_exon_intron_ratios.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/QC_hg38/K562_RNA-seq_90_exon_intron_ratios.bed` (110022)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.53GB.  
  PID: 110022;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Produce bigWig files (02-18 11:39:53) elapsed: 41.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/signal_hg38/K562_RNA-seq_90_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/signal_hg38/K562_RNA-seq_90_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_plus.bam` (110029)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.53GB.  
  PID: 110029;	Command: samtools;	Return code: 0;	Memory used: 0.005GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/signal_hg38/K562_RNA-seq_90_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/signal_hg38/K562_RNA-seq_90_plus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge` (110033)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_90_plus_cuttrace_anzaahji'
Processing with 1 cores...
Discarding 131 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrEBV']
Keeping 64 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270750v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 64 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/signal_hg38/K562_RNA-seq_90_plus_exact_body_0-mer.bw'
Merging 64 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/signal_hg38/K562_RNA-seq_90_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:14:47. Running peak memory: 3.53GB.  
  PID: 110033;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.885GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/signal_hg38/K562_RNA-seq_90_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/signal_hg38/K562_RNA-seq_90_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_minus.bam` (113361)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.53GB.  
  PID: 113361;	Command: samtools;	Return code: 0;	Memory used: 0.005GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/signal_hg38/K562_RNA-seq_90_minus_exact_body_0-mer.bw -p 2 --variable-step --tail-edge` (113364)
<pre>
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/aligned_hg38/K562_RNA-seq_90_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_90_minus_cuttrace_cx9u6q12'
Processing with 2 cores...
Discarding 133 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr16_KI270728v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 62 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr17_GL000205v2_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 62 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_RNA-seq_90/signal_hg38/K562_RNA-seq_90_minus_exact_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.53GB.  
  PID: 113364;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.102GB

Starting cleanup: 52 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:46:01
*  Total elapsed time (all runs):  1:04:36
*         Peak memory (this run):  3.5297 GB
*        Pipeline completed time: 2020-02-18 11:54:58
