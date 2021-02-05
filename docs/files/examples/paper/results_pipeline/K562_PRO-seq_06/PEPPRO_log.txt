### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_06 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_6pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty`
*         Compute host:  udc-ba26-34c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/
*  Pipeline started at:   (06-11 16:59:53) elapsed: 7.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_6pct.fastq.gz']`
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
*        `sample_name`:  `K562_PRO-seq_06`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_6pct.fastq.gz

> `File_mb`	2220.82	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 16:59:54) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/K562_PRO-seq_06.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_PRO_6pct.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/K562_PRO-seq_06.fastq.gz` (170884)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 170884;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/K562_PRO-seq_06.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/K562_PRO-seq_06_R1.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/K562_PRO-seq_06.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/K562_PRO-seq_06_R1.fastq` (170885)
<pre>
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 0.002GB.  
  PID: 170885;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	29786346	PEPPRO	_RES_

> `Fastq_reads`	29786346	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/K562_PRO-seq_06.fastq.gz']

### FASTQ processing:  (06-11 17:01:25) elapsed: 91.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/K562_PRO-seq_06_R1_processed.fastq`  

> `(cutadapt -j 12 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/K562_PRO-seq_06_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/K562_PRO-seq_06_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/K562_PRO-seq_06_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt/K562_PRO-seq_06_R1_cutadapt.txt` (171343)
<pre>
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 3.955GB.  
  PID: 171343;	Command: cutadapt;	Return code: 0;	Memory used: 3.955GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/K562_PRO-seq_06_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/K562_PRO-seq_06_R1_processed.fastq` (171415,171416)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 3.955GB.  
  PID: 171415;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 171416;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt/K562_PRO-seq_06_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	25373082.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt/K562_PRO-seq_06_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt/K562_PRO-seq_06_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/K562_PRO-seq_06_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	747032.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	2.508	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastqc/K562_PRO-seq_06_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (171938)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.955GB.  
  PID: 171938;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	29039314	PEPPRO	_RES_

> `Trim_loss_rate`	2.51	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/K562_PRO-seq_06_R1_processed.fastq` (171956)
<pre>
Started analysis of K562_PRO-seq_06_R1_processed.fastq
Approx 5% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 10% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 15% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 20% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 25% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 30% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 35% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 40% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 45% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 50% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 55% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 60% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 65% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 70% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 75% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 80% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 85% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 90% complete for K562_PRO-seq_06_R1_processed.fastq
Approx 95% complete for K562_PRO-seq_06_R1_processed.fastq
Analysis complete for K562_PRO-seq_06_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:40. Running peak memory: 3.955GB.  
  PID: 171956;	Command: fastqc;	Return code: 0;	Memory used: 0.197GB

> `FastQC report r1`	fastqc/K562_PRO-seq_06_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/processed_R1.flag` (172162)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.955GB.  
  PID: 172162;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-11 17:04:48) elapsed: 203.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt/K562_PRO-seq_06_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt/K562_PRO-seq_06_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt` (172163)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.955GB.  
  PID: 172163;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_06_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_06_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt/K562_PRO-seq_06_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 17:04:54) elapsed: 6.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt/K562_PRO-seq_06_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt/K562_PRO-seq_06_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt/K562_PRO-seq_06_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt/K562_PRO-seq_06_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt/K562_PRO-seq_06_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2311	PEPPRO	_RES_

### Prealignments (06-11 17:04:54) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 17:04:54) elapsed: 0.0 _TIME_


> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_PRO-seq_06 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/K562_PRO-seq_06_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/prealignments/K562_PRO-seq_06_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
29039314 reads; of these:
  29039314 (100.00%) were unpaired; of these:
    26359027 (90.77%) aligned 0 times
    2680287 (9.23%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.23% overall alignment rate

> `Aligned_reads_human_rDNA`	2680287.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.23	PEPPRO	_RES_

### Map to genome (06-11 17:07:54) elapsed: 180.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam`  

> `bowtie2 -p 12 --very-sensitive --rg-id K562_PRO-seq_06 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/prealignments/K562_PRO-seq_06_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/tmpkd1ui_rj -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_temp.bam` (172819,172820,172821)
<pre>
26359027 reads; of these:
  26359027 (100.00%) were unpaired; of these:
    324563 (1.23%) aligned 0 times
    19683305 (74.67%) aligned exactly 1 time
    6351159 (24.09%) aligned >1 times
98.77% overall alignment rate
[bam_sort_core] merging from 8 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:18:56. Running peak memory: 3.955GB.  
  PID: 172819;	Command: bowtie2;	Return code: 0;	Memory used: 3.717GB  
  PID: 172820;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 172821;	Command: samtools;	Return code: 0;	Memory used: 0.888GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam` (178942)
<pre>
</pre>
Command completed. Elapsed time: 0:00:56. Running peak memory: 3.955GB.  
  PID: 178942;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	26034464	PEPPRO	_RES_

> `QC_filtered_reads`	2801424	PEPPRO	_RES_

> `Aligned_reads`	23233040	PEPPRO	_RES_

> `Alignment_rate`	80.01	PEPPRO	_RES_

> `Total_efficiency`	78.0	PEPPRO	_RES_

> `Read_depth`	3.34	PEPPRO	_RES_

### Compress all unmapped read files (06-11 17:42:28) elapsed: 2074.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/prealignments/K562_PRO-seq_06_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/prealignments/K562_PRO-seq_06_human_rDNA_unmap.fq` (182120)
<pre>
</pre>
Command completed. Elapsed time: 0:00:50. Running peak memory: 3.955GB.  
  PID: 182120;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_temp.bam` (182193)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 3.955GB.  
  PID: 182193;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	548203	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam` (182222)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 3.955GB.  
  PID: 182222;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/chr_sizes.bed` (182255,182256,182257,182258)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.955GB.  
  PID: 182255;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 182257;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 182256;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 182258;	Command: grep;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_noMT.bam` (182260)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.955GB.  
  PID: 182260;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam` (182508)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.955GB.  
  PID: 182508;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam` (182510)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 3.955GB.  
  PID: 182510;	Command: samtools;	Return code: 0;	Memory used: 0.011GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-11 17:46:21) elapsed: 233.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_bamQC.tsv` (182754)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/tmp_K562_PRO-seq_06_sort_5arlm633'
Processing with 12 cores...
Discarding 100 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 95 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 3.955GB.  
  PID: 182754;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.496GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_bamQC.tsv`

> `NRF`	0.87	PEPPRO	_RES_

> `PBC1`	0.94	PEPPRO	_RES_

> `PBC2`	19.95	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_unmap.bam`  

> `samtools view -b -@ 12 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_unmap.bam` (182825)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.955GB.  
  PID: 182825;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_temp.bam`

> `Unmapped_reads`	324563	PEPPRO	_RES_

### Split BAM by strand (06-11 17:46:58) elapsed: 37.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam` (182861)
<pre>
</pre>
Command completed. Elapsed time: 0:01:53. Running peak memory: 3.955GB.  
  PID: 182861;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam` (182975)
<pre>
</pre>
Command completed. Elapsed time: 0:01:50. Running peak memory: 3.955GB.  
  PID: 182975;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-11 17:50:41) elapsed: 224.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (183542)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.955GB.  
  PID: 183542;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_plus_TssEnrichment.txt` (183544)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.955GB.  
  PID: 183544;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.709GB


> `TSS_coding_score`	13.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_minus_TssEnrichment.txt` (183575)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.955GB.  
  PID: 183575;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.814GB


> `TSS_non-coding_score`	5.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_minus_TssEnrichment.txt` (183606)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.955GB.  
  PID: 183606;	Command: Rscript;	Return code: 0;	Memory used: 0.222GB

> `TSS enrichment`	QC_hg38/K562_PRO-seq_06_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_06_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt` (183631,183632,183633,183634)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.955GB.  
  PID: 183631;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 183633;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 183632;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 183634;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_keep.txt` (183636)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.955GB.  
  PID: 183636;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-11 17:51:01) elapsed: 19.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_ensembl_tss.bed` (183639,183640)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.955GB.  
  PID: 183639;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 183640;	Command: bedtools;	Return code: 0;	Memory used: 0.1GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_ensembl_gene_body.bed` (183643,183644)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.955GB.  
  PID: 183643;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 183644;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_TSS_density.bed` (183646,183647,183648,183649)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.955GB.  
  PID: 183647;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 183649;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 183646;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 183648;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_gene_body_density.bed` (183682,183683,183684)
<pre>
</pre>
Command completed. Elapsed time: 0:00:43. Running peak memory: 3.955GB.  
  PID: 183682;	Command: bedtools;	Return code: 0;	Memory used: 0.084GB  
  PID: 183684;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 183683;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_pause_index.bed`  
awk: 
> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/tmp67q5ndox` (183721,183722,183723)
cmd. line:1: { if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}
awk: cmd. line:1:             ^ syntax error
awk: cmd. line:1: { if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}
awk: cmd. line:1:                                                                                                                      ^ syntax error
<pre>
</pre>
Got SIGTERM. Failing gracefully... (06-11 18:55:13) elapsed: 3853.0 _TIME_
Child process 183721 (join) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/recover.lock.QC_hg38__K562_PRO-seq_06_pause_index.bed

### Pipeline failed at:  (06-11 18:55:13) elapsed: 0.0 _TIME_

Total time: 1:55:27
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_06 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_6pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aw29-25a
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/
*  Pipeline started at:   (06-11 19:03:04) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_6pct.fastq.gz']`
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
*        `sample_name`:  `K562_PRO-seq_06`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_6pct.fastq.gz

> `File_mb`	2220.82	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 19:03:05) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/K562_PRO-seq_06.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/K562_PRO-seq_06.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/K562_PRO-seq_06_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/K562_PRO-seq_06.fastq.gz']

### FASTQ processing:  (06-11 19:03:05) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-11 19:03:05) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt/K562_PRO-seq_06_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_06_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_06_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-11 19:03:05) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 19:03:05) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/prealignments/K562_PRO-seq_06_human_rDNA_unmap.fq

### Map to genome (06-11 19:03:05) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam`  

### Compress all unmapped read files (06-11 19:03:05) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/prealignments/K562_PRO-seq_06_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-11 19:03:05) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_unmap.bam`  

### Split BAM by strand (06-11 19:03:05) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam`  

### Calculate TSS enrichment (06-11 19:03:05) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/K562_PRO-seq_06_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_06_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-11 19:03:05) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_gene_body_density.bed`  
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/lock.QC_hg38__K562_PRO-seq_06_pause_index.bed
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/tmpddr5l236` (72731,72732,72733)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 72731;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 72733;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 72732;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/tmpddr5l236 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-11 19:03:05) elapsed: 0.0 _TIME_

Total time: 0:00:01
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/tmpddr5l236 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_06 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_6pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba26-33c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/
*  Pipeline started at:   (06-14 21:08:16) elapsed: 124.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_6pct.fastq.gz']`
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
*        `sample_name`:  `K562_PRO-seq_06`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_6pct.fastq.gz

> `File_mb`	2220.82	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-14 21:08:17) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/K562_PRO-seq_06.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/K562_PRO-seq_06.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/K562_PRO-seq_06_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/K562_PRO-seq_06.fastq.gz']

### FASTQ processing:  (06-14 21:08:17) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-14 21:08:17) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/cutadapt/K562_PRO-seq_06_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_06_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_06_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-14 21:08:17) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:08:17) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/prealignments/K562_PRO-seq_06_human_rDNA_unmap.fq

### Map to genome (06-14 21:08:17) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam`  

### Compress all unmapped read files (06-14 21:08:17) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/prealignments/K562_PRO-seq_06_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:08:17) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_unmap.bam`  

### Split BAM by strand (06-14 21:08:17) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam`  

### Calculate TSS enrichment (06-14 21:08:17) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/K562_PRO-seq_06_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_06_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:08:17) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/tmplens6kmt` (1547,1548,1549)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 1547;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 1549;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 1548;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/tmplens6kmt | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.02107) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/tmplens6kmt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_pause_index.bed` (1556)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 1556;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.67	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_pause_index.bed` (1561)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.204GB.  
  PID: 1561;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Pause index`	QC_hg38/K562_PRO-seq_06_pause_index.pdf	Pause index	QC_hg38/K562_PRO-seq_06_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_pause_index.bed` (1587)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 1587;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:08:24) elapsed: 7.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam`
23233040 8290569

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam`
23233040 7882357

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/signal_hg38/K562_PRO-seq_06_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_gene_sort.bed` (1640,1641)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 1640;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 1641;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/signal_hg38/K562_PRO-seq_06_gene_coverage.bed` (1644)
<pre>
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 0.204GB.  
  PID: 1644;	Command: bedtools;	Return code: 0;	Memory used: 0.082GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/hg38_annotations.bed.gz` (1689)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 1689;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/hg38_annotations.bed` (1690)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 1690;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:09:57) elapsed: 92.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/raw/hg38_annotations.bed` (1699)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 1699;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Enhancer_sort.bed` (1702,1703,1704,1705)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 1702;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 1703;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 1705;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB  
  PID: 1704;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Enhancer_plus_coverage.bed` (1707)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 0.204GB.  
  PID: 1707;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Enhancer_minus_coverage.bed` (1971)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 0.204GB.  
  PID: 1971;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter_sort.bed` (1990,1991,1992,1993)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 1990;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 1991;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 1993;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 1992;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Promoter_plus_coverage.bed` (1995)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 0.204GB.  
  PID: 1995;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Promoter_minus_coverage.bed` (2014)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 0.204GB.  
  PID: 2014;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter_Flanking_Region"` (2032)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 2032;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter_Flanking_Region_sort.bed` (2034,2035,2036,2037)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 2034;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 2036;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 2035;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 2037;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Promoter_Flanking_Region_plus_coverage.bed` (2039)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 0.204GB.  
  PID: 2039;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Promoter_Flanking_Region_minus_coverage.bed` (2060)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 0.204GB.  
  PID: 2060;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/5_UTR"` (2077)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 2077;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/5_UTR_sort.bed` (2078,2079,2080,2081)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 2078;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 2079;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 2081;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 2080;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_5_UTR_plus_coverage.bed` (2083)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 0.204GB.  
  PID: 2083;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_5_UTR_minus_coverage.bed` (2100)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 0.204GB.  
  PID: 2100;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/3_UTR"` (2119)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 2119;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/3_UTR_sort.bed` (2120,2121,2122,2123)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 2120;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 2121;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 2123;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 2122;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_3_UTR_plus_coverage.bed` (2125)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 0.204GB.  
  PID: 2125;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_3_UTR_minus_coverage.bed` (2141)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 0.204GB.  
  PID: 2141;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Exon_sort.bed` (2157,2158,2159,2160)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 2157;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 2158;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 2160;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 2159;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Exon_plus_coverage.bed` (2165)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 0.204GB.  
  PID: 2165;	Command: bedtools;	Return code: 0;	Memory used: 0.038GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Exon_minus_coverage.bed` (2186)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 0.204GB.  
  PID: 2186;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Intron_sort.bed` (2203,2205,2206,2207)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 2203;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 2206;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 2205;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 2207;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Intron_plus_coverage.bed` (2210)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 0.204GB.  
  PID: 2210;	Command: bedtools;	Return code: 0;	Memory used: 0.029GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Intron_minus_coverage.bed` (2230)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 0.204GB.  
  PID: 2230;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB


### Plot cFRiF/FRiF (06-14 21:14:28) elapsed: 272.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_06 -z 3099922541 -n 11562878 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Intron_plus_coverage.bed` (2268)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 0.456GB.  
  PID: 2268;	Command: Rscript;	Return code: 0;	Memory used: 0.456GB

> `cFRiF`	QC_hg38/K562_PRO-seq_06_cFRiF.pdf	cFRiF	QC_hg38/K562_PRO-seq_06_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_06 -z 3099922541 -n 11562878 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_Intron_plus_coverage.bed` (2501)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 0.456GB.  
  PID: 2501;	Command: Rscript;	Return code: 0;	Memory used: 0.456GB

> `FRiF`	QC_hg38/K562_PRO-seq_06_FRiF.pdf	FRiF	QC_hg38/K562_PRO-seq_06_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 21:15:34) elapsed: 66.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_exons_sort.bed` (2542,2543)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.456GB.  
  PID: 2543;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB  
  PID: 2542;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_introns_sort.bed` (2551,2552,2553)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.456GB.  
  PID: 2551;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 2553;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 2552;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_exons_coverage.bed` (2560)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 0.456GB.  
  PID: 2560;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_introns_coverage.bed` (2595)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 0.456GB.  
  PID: 2595;	Command: bedtools;	Return code: 0;	Memory used: 0.032GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/23.23304)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_exons_rpkm.bed` (2629,2630,2631)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.456GB.  
  PID: 2629;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 2631;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 2630;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/23.23304)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_introns_rpkm.bed` (2633,2634,2635)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.456GB.  
  PID: 2633;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 2635;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 2634;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_exon_intron_ratios.bed` (2638,2639,2640)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.456GB.  
  PID: 2638;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 2640;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 2639;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.41	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_exon_intron_ratios.bed --annotate` (2646)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.456GB.  
  PID: 2646;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/K562_PRO-seq_06_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_PRO-seq_06_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/QC_hg38/K562_PRO-seq_06_exon_intron_ratios.bed` (2667)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.456GB.  
  PID: 2667;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-14 21:17:06) elapsed: 92.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/signal_hg38/K562_PRO-seq_06_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/signal_hg38/K562_PRO-seq_06_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam` (2675)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 0.456GB.  
  PID: 2675;	Command: samtools;	Return code: 0;	Memory used: 0.007GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/signal_hg38/K562_PRO-seq_06_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/signal_hg38/K562_PRO-seq_06_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 23233040.0` (2685)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_plus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_06_plus_cuttrace_cuoiph9b'
Processing with 4 cores...
Discarding 114 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr22_KI270732v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrEBV']
Keeping 81 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 81 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/signal_hg38/K562_PRO-seq_06_plus_exact_body_0-mer.bw'
Merging 81 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/signal_hg38/K562_PRO-seq_06_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:09. Running peak memory: 2.621GB.  
  PID: 2685;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.621GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/signal_hg38/K562_PRO-seq_06_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/signal_hg38/K562_PRO-seq_06_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam` (4340)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 2.621GB.  
  PID: 4340;	Command: samtools;	Return code: 0;	Memory used: 0.007GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/signal_hg38/K562_PRO-seq_06_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/signal_hg38/K562_PRO-seq_06_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 23233040.0` (4355)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/aligned_hg38/K562_PRO-seq_06_minus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_06_minus_cuttrace_g8q_7pyu'
Processing with 4 cores...
stdin is empty of data
Discarding 112 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1']
Keeping 83 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/signal_hg38/K562_PRO-seq_06_minus_exact_body_0-mer.bw'
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_06/signal_hg38/K562_PRO-seq_06_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:53. Running peak memory: 2.645GB.  
  PID: 4355;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.645GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:27:19
*  Total elapsed time (all runs):  1:27:32
*         Peak memory (this run):  2.6445 GB
*        Pipeline completed time: 2020-06-14 21:33:32
