### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_30 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_30pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty`
*         Compute host:  udc-aj38-16c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/
*  Pipeline started at:   (06-11 17:09:22) elapsed: 2.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_30pctRNA.fastq.gz']`
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
*            `recover`:  `False`
*        `sample_name`:  `K562_RNA-seq_30`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_30pctRNA.fastq.gz

> `File_mb`	5620.96	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 17:09:30) elapsed: 8.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_30pctRNA.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz` (455082)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 455082;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1.fastq` (455084)
<pre>
</pre>
Command completed. Elapsed time: 0:03:24. Running peak memory: 0.002GB.  
  PID: 455084;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	70000000	PEPPRO	_RES_

> `Fastq_reads`	70000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz']

### FASTQ processing:  (06-11 17:18:32) elapsed: 542.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_processed.fastq`  

> `(cutadapt -j 12 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt` (456083)
<pre>
</pre>
Command completed. Elapsed time: 0:04:32. Running peak memory: 4.222GB.  
  PID: 456083;	Command: cutadapt;	Return code: 0;	Memory used: 4.222GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_processed.fastq` (456630,456631)
<pre>
</pre>
Command completed. Elapsed time: 0:01:46. Running peak memory: 4.222GB.  
  PID: 456630;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 456631;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	49200456.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	1225550.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	1.7508	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastqc/K562_RNA-seq_30_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (457078)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.222GB.  
  PID: 457078;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	68774450	PEPPRO	_RES_

> `Trim_loss_rate`	1.75	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_processed.fastq` (457146)
<pre>
Started analysis of K562_RNA-seq_30_R1_processed.fastq
Approx 5% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 10% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 15% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 20% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 25% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 30% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 35% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 40% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 45% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 50% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 55% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 60% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 65% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 70% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 75% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 80% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 85% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 90% complete for K562_RNA-seq_30_R1_processed.fastq
Approx 95% complete for K562_RNA-seq_30_R1_processed.fastq
Analysis complete for K562_RNA-seq_30_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:04:06. Running peak memory: 4.222GB.  
  PID: 457146;	Command: fastqc;	Return code: 0;	Memory used: 0.242GB

> `FastQC report r1`	fastqc/K562_RNA-seq_30_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/processed_R1.flag` (457691)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.222GB.  
  PID: 457691;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Plot adapter insertion distribution (06-11 17:32:08) elapsed: 816.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt` (457692)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 4.222GB.  
  PID: 457692;	Command: Rscript;	Return code: 0;	Memory used: 0.201GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 17:32:17) elapsed: 9.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.231	PEPPRO	_RES_

### Prealignments (06-11 17:32:17) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 17:32:17) elapsed: 0.0 _TIME_


> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_30 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
68774450 reads; of these:
  68774450 (100.00%) were unpaired; of these:
    63879577 (92.88%) aligned 0 times
    4894873 (7.12%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
7.12% overall alignment rate

> `Aligned_reads_human_rDNA`	4894873.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	7.12	PEPPRO	_RES_

### Map to genome (06-11 17:41:27) elapsed: 550.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam`  

> `bowtie2 -p 12 --very-sensitive --rg-id K562_RNA-seq_30 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/tmpt4mjstb1 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam` (300,308,309)
<pre>
63879577 reads; of these:
  63879577 (100.00%) were unpaired; of these:
    3991241 (6.25%) aligned 0 times
    42466485 (66.48%) aligned exactly 1 time
    17421851 (27.27%) aligned >1 times
93.75% overall alignment rate
[bam_sort_core] merging from 20 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:41:45. Running peak memory: 4.222GB.  
  PID: 308;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 300;	Command: bowtie2;	Return code: 0;	Memory used: 3.718GB  
  PID: 309;	Command: samtools;	Return code: 0;	Memory used: 0.889GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam` (10971)
<pre>
</pre>
Command completed. Elapsed time: 0:02:34. Running peak memory: 4.222GB.  
  PID: 10971;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	59888336	PEPPRO	_RES_

> `QC_filtered_reads`	8412839	PEPPRO	_RES_

> `Aligned_reads`	51475497	PEPPRO	_RES_

> `Alignment_rate`	74.85	PEPPRO	_RES_

> `Total_efficiency`	73.54	PEPPRO	_RES_

> `Read_depth`	6.15	PEPPRO	_RES_

### Compress all unmapped read files (06-11 18:46:27) elapsed: 3900.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq` (19308)
<pre>
</pre>
Command completed. Elapsed time: 0:02:29. Running peak memory: 4.222GB.  
  PID: 19308;	Command: pigz;	Return code: 0;	Memory used: 0.014GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam` (20001)
<pre>
</pre>
Command completed. Elapsed time: 0:01:04. Running peak memory: 4.222GB.  
  PID: 20001;	Command: samtools;	Return code: 0;	Memory used: 0.017GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	2160789	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam` (20325)
<pre>
</pre>
Command completed. Elapsed time: 0:00:52. Running peak memory: 4.222GB.  
  PID: 20325;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/chr_sizes.bed` (20782,20783,20784,20785)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.222GB.  
  PID: 20783;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 20785;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 20782;	Command: samtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 20784;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_noMT.bam` (20787)
<pre>
</pre>
Command completed. Elapsed time: 0:01:04. Running peak memory: 4.222GB.  
  PID: 20787;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam` (20922)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.222GB.  
  PID: 20922;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam` (20923)
<pre>
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 4.222GB.  
  PID: 20923;	Command: samtools;	Return code: 0;	Memory used: 0.015GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`
</pre>
Got SIGTERM. Failing gracefully... (06-11 18:55:14) elapsed: 527.0 _TIME_
No locked process. Dynamic recovery will be automatic.

### Pipeline failed at:  (06-11 18:55:14) elapsed: 0.0 _TIME_

Total time: 1:45:55
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_30 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_30pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aw29-25b
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/
*  Pipeline started at:   (06-11 19:08:12) elapsed: 2.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_30pctRNA.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_30`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_30pctRNA.fastq.gz

> `File_mb`	5620.96	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 19:08:13) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz']

### FASTQ processing:  (06-11 19:08:13) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-11 19:08:13) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-11 19:08:13) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 19:08:13) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq

### Map to genome (06-11 19:08:13) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam`  

### Compress all unmapped read files (06-11 19:08:13) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq.gz`  
Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-11 19:10:57) elapsed: 164.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv` (414538)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/tmp_K562_RNA-seq_30_sort_yz6hm5_8'
Processing with 12 cores...
Discarding 97 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 98 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:01:10. Running peak memory: 1.422GB.  
  PID: 414538;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.422GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv`

> `NRF`	0.77	PEPPRO	_RES_

> `PBC1`	0.9	PEPPRO	_RES_

> `PBC2`	13.32	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_unmap.bam`  

> `samtools view -b -@ 12 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_unmap.bam` (414947)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 1.422GB.  
  PID: 414947;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam`

> `Unmapped_reads`	3991241	PEPPRO	_RES_

### Split BAM by strand (06-11 19:12:35) elapsed: 98.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam` (415009)
<pre>
</pre>
Command completed. Elapsed time: 0:04:06. Running peak memory: 1.422GB.  
  PID: 415009;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam` (415553)
<pre>
</pre>
Command completed. Elapsed time: 0:03:57. Running peak memory: 1.422GB.  
  PID: 415553;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-11 19:20:39) elapsed: 484.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (416344)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.422GB.  
  PID: 416344;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_plus_TssEnrichment.txt` (416347)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 1.422GB.  
  PID: 416347;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.533GB


> `TSS_coding_score`	16.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_minus_TssEnrichment.txt` (416559)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 1.422GB.  
  PID: 416559;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.543GB


> `TSS_non-coding_score`	4.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_minus_TssEnrichment.txt` (416594)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 1.422GB.  
  PID: 416594;	Command: Rscript;	Return code: 0;	Memory used: 0.206GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_30_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_30_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt` (416619,416620,416621,416622)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.422GB.  
  PID: 416619;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 416621;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 416620;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 416622;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt` (416625)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.422GB.  
  PID: 416625;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-11 19:21:05) elapsed: 27.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_tss.bed` (416627,416628)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 1.422GB.  
  PID: 416627;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 416628;	Command: bedtools;	Return code: 0;	Memory used: 0.096GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_gene_body.bed` (416631,416632)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.422GB.  
  PID: 416631;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 416632;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSS_density.bed` (416634,416635,416636,416637)
<pre>
</pre>
Command completed. Elapsed time: 0:01:17. Running peak memory: 1.422GB.  
  PID: 416635;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 416637;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 416634;	Command: bedtools;	Return code: 0;	Memory used: 0.026GB  
  PID: 416636;	Command: sort;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_gene_body_density.bed` (416706,416707,416708)
<pre>
</pre>
Command completed. Elapsed time: 0:01:52. Running peak memory: 1.422GB.  
  PID: 416706;	Command: bedtools;	Return code: 0;	Memory used: 0.15GB  
  PID: 416708;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 416707;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/tmp5xxw8i09` (416807,416808,416809)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.422GB.  
  PID: 416807;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 416809;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 416808;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/tmp5xxw8i09 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-11 19:24:17) elapsed: 192.0 _TIME_

Total time: 0:16:07
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/tmp5xxw8i09 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_30 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_30pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba25-13c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/
*  Pipeline started at:   (06-14 21:11:24) elapsed: 7.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_30pctRNA.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_30`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_30pctRNA.fastq.gz

> `File_mb`	5620.96	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-14 21:11:25) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz']

### FASTQ processing:  (06-14 21:11:25) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-14 21:11:25) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-14 21:11:25) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:11:25) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq

### Map to genome (06-14 21:11:25) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam`  

### Compress all unmapped read files (06-14 21:11:25) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:11:25) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_unmap.bam`  

### Split BAM by strand (06-14 21:11:25) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam`  

### Calculate TSS enrichment (06-14 21:11:25) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/K562_RNA-seq_30_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_30_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:11:25) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/tmppapopl94` (409549,409550,409551)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 409549;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 409551;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 409550;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/tmppapopl94 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0417111) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/tmppapopl94 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed` (409557)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 409557;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	8.95	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed` (409562)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.316GB.  
  PID: 409562;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Pause index`	QC_hg38/K562_RNA-seq_30_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_30_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed` (409583)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.316GB.  
  PID: 409583;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:11:30) elapsed: 5.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam`
51475497 18930658

> `Plus_FRiP`	0.37	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam`
51475497 18369003

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_gene_sort.bed` (409678,409679)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.316GB.  
  PID: 409678;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 409679;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_gene_coverage.bed` (409682)
<pre>
</pre>
Command completed. Elapsed time: 0:02:19. Running peak memory: 0.316GB.  
  PID: 409682;	Command: bedtools;	Return code: 0;	Memory used: 0.137GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed.gz` (410001)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.316GB.  
  PID: 410001;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed` (410002)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.316GB.  
  PID: 410002;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:15:20) elapsed: 230.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed` (410013)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.316GB.  
  PID: 410013;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer_sort.bed` (410015,410016,410017,410018)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.316GB.  
  PID: 410015;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 410016;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 410018;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB  
  PID: 410017;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_plus_coverage.bed` (410020)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 0.316GB.  
  PID: 410020;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_minus_coverage.bed` (410060)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 0.316GB.  
  PID: 410060;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_sort.bed` (410095,410096,410097,410098)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.316GB.  
  PID: 410095;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 410096;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 410098;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 410097;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_plus_coverage.bed` (410100)
<pre>
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 0.316GB.  
  PID: 410100;	Command: bedtools;	Return code: 0;	Memory used: 0.049GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_minus_coverage.bed` (410140)
<pre>
</pre>
Command completed. Elapsed time: 0:00:44. Running peak memory: 0.316GB.  
  PID: 410140;	Command: bedtools;	Return code: 0;	Memory used: 0.062GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region"` (410180)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.316GB.  
  PID: 410180;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed` (410181,410182,410183,410184)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.316GB.  
  PID: 410181;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 410183;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 410182;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 410184;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_plus_coverage.bed` (410187)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 0.316GB.  
  PID: 410187;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_minus_coverage.bed` (410223)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 0.316GB.  
  PID: 410223;	Command: bedtools;	Return code: 0;	Memory used: 0.027GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR"` (410260)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.316GB.  
  PID: 410260;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR_sort.bed` (410261,410262,410263,410264)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.316GB.  
  PID: 410261;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 410262;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 410264;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 410263;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_plus_coverage.bed` (410267)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 0.316GB.  
  PID: 410267;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_minus_coverage.bed` (410519)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 0.316GB.  
  PID: 410519;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR"` (410576)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.316GB.  
  PID: 410576;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR_sort.bed` (410578,410579,410580,410581)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.316GB.  
  PID: 410578;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 410579;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 410581;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB  
  PID: 410580;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_plus_coverage.bed` (410583)
<pre>
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 0.316GB.  
  PID: 410583;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_minus_coverage.bed` (410623)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 0.316GB.  
  PID: 410623;	Command: bedtools;	Return code: 0;	Memory used: 0.032GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon_sort.bed` (410660,410661,410662,410663)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.316GB.  
  PID: 410660;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 410661;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 410663;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 410662;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_plus_coverage.bed` (410667)
<pre>
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 0.316GB.  
  PID: 410667;	Command: bedtools;	Return code: 0;	Memory used: 0.091GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_minus_coverage.bed` (410707)
<pre>
</pre>
Command completed. Elapsed time: 0:00:44. Running peak memory: 0.316GB.  
  PID: 410707;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron_sort.bed` (410747,410748,410749,410750)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.316GB.  
  PID: 410747;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 410749;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 410748;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 410750;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_plus_coverage.bed` (410754)
<pre>
</pre>
Command completed. Elapsed time: 0:00:48. Running peak memory: 0.316GB.  
  PID: 410754;	Command: bedtools;	Return code: 0;	Memory used: 0.075GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_minus_coverage.bed` (410797)
<pre>
</pre>
Command completed. Elapsed time: 0:00:50. Running peak memory: 0.316GB.  
  PID: 410797;	Command: bedtools;	Return code: 0;	Memory used: 0.063GB


### Plot cFRiF/FRiF (06-14 21:25:18) elapsed: 598.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_RNA-seq_30 -z 3099922541 -n 25105669 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_plus_coverage.bed` (411051)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 0.456GB.  
  PID: 411051;	Command: Rscript;	Return code: 0;	Memory used: 0.456GB

> `cFRiF`	QC_hg38/K562_RNA-seq_30_cFRiF.pdf	cFRiF	QC_hg38/K562_RNA-seq_30_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_RNA-seq_30 -z 3099922541 -n 25105669 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_plus_coverage.bed` (411104)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 0.458GB.  
  PID: 411104;	Command: Rscript;	Return code: 0;	Memory used: 0.458GB

> `FRiF`	QC_hg38/K562_RNA-seq_30_FRiF.pdf	FRiF	QC_hg38/K562_RNA-seq_30_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 21:26:34) elapsed: 76.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_exons_sort.bed` (411139,411140)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.458GB.  
  PID: 411140;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB  
  PID: 411139;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_introns_sort.bed` (411147,411148,411149)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.458GB.  
  PID: 411147;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 411149;	Command: bedtools;	Return code: 0;	Memory used: 0.003GB  
  PID: 411148;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_coverage.bed` (411156)
<pre>
</pre>
Command completed. Elapsed time: 0:01:28. Running peak memory: 0.458GB.  
  PID: 411156;	Command: bedtools;	Return code: 0;	Memory used: 0.074GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_coverage.bed` (411235)
<pre>
</pre>
Command completed. Elapsed time: 0:01:37. Running peak memory: 0.458GB.  
  PID: 411235;	Command: bedtools;	Return code: 0;	Memory used: 0.07GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/51.475497)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_rpkm.bed` (411321,411322,411323)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.458GB.  
  PID: 411321;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 411323;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 411322;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/51.475497)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_rpkm.bed` (411326,411327,411328)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.458GB.  
  PID: 411326;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 411328;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 411327;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed` (411330,411331,411332)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.458GB.  
  PID: 411330;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 411332;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 411331;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	2.97	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed --annotate` (411339)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.458GB.  
  PID: 411339;	Command: Rscript;	Return code: 0;	Memory used: 0.285GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_30_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_30_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed` (411360)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.458GB.  
  PID: 411360;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-14 21:29:57) elapsed: 203.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam` (411369)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 0.458GB.  
  PID: 411369;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 51475497.0` (411649)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_30_plus_cuttrace_ou77znyr'
Processing with 4 cores...
stdin is empty of data
Discarding 109 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 86 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_exact_body_0-mer.bw'
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:33. Running peak memory: 3.453GB.  
  PID: 411649;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.453GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam` (413141)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 3.453GB.  
  PID: 413141;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 51475497.0` (413165)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_30_minus_cuttrace_f3x6q4jg'
Processing with 4 cores...
stdin is empty of data
Discarding 110 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1']
Keeping 85 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 85 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_exact_body_0-mer.bw'
Merging 85 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:26. Running peak memory: 3.453GB.  
  PID: 413165;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.244GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:36:29
*  Total elapsed time (all runs):  2:47:25
*         Peak memory (this run):  3.453 GB
*        Pipeline completed time: 2020-06-14 21:47:46
