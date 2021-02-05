### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_10 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_10pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty`
*         Compute host:  udc-ba27-10
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/
*  Pipeline started at:   (06-11 17:08:25) elapsed: 3.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_10pctRNA.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_10`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_10pctRNA.fastq.gz

> `File_mb`	5660.06	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 17:08:26) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_10pctRNA.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz` (58545)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 58545;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1.fastq` (58546)
<pre>
</pre>
Command completed. Elapsed time: 0:03:02. Running peak memory: 0.002GB.  
  PID: 58546;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	70000000	PEPPRO	_RES_

> `Fastq_reads`	70000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz']

### FASTQ processing:  (06-11 17:13:07) elapsed: 280.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_processed.fastq`  

> `(cutadapt -j 12 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt` (59233)
<pre>
</pre>
Command completed. Elapsed time: 0:02:26. Running peak memory: 4.589GB.  
  PID: 59233;	Command: cutadapt;	Return code: 0;	Memory used: 4.589GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_processed.fastq` (59614,59615)
<pre>
</pre>
Command completed. Elapsed time: 0:02:23. Running peak memory: 4.589GB.  
  PID: 59615;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB  
  PID: 59614;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	56155013.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	1576345.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	2.2519	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastqc/K562_RNA-seq_10_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (61043)
### Calculate the number of trimmed reads
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.589GB.  
  PID: 61043;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	68423655	PEPPRO	_RES_

> `Trim_loss_rate`	2.25	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_processed.fastq` (61201)
<pre>
Started analysis of K562_RNA-seq_10_R1_processed.fastq
Approx 5% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 10% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 15% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 20% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 25% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 30% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 35% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 40% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 45% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 50% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 55% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 60% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 65% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 70% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 75% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 80% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 85% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 90% complete for K562_RNA-seq_10_R1_processed.fastq
Approx 95% complete for K562_RNA-seq_10_R1_processed.fastq
Analysis complete for K562_RNA-seq_10_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:04:06. Running peak memory: 4.589GB.  
  PID: 61201;	Command: fastqc;	Return code: 0;	Memory used: 0.238GB

> `FastQC report r1`	fastqc/K562_RNA-seq_10_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/processed_R1.flag` (63276)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.589GB.  
  PID: 63276;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-11 17:23:50) elapsed: 643.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt` (63277)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 4.589GB.  
  PID: 63277;	Command: Rscript;	Return code: 0;	Memory used: 0.111GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 17:23:59) elapsed: 9.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.231	PEPPRO	_RES_

### Prealignments (06-11 17:23:59) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 17:23:59) elapsed: 0.0 _TIME_


> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_10 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
68423655 reads; of these:
  68423655 (100.00%) were unpaired; of these:
    62592060 (91.48%) aligned 0 times
    5831595 (8.52%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
8.52% overall alignment rate

> `Aligned_reads_human_rDNA`	5831595.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	8.52	PEPPRO	_RES_

### Map to genome (06-11 17:31:41) elapsed: 462.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam`  

> `bowtie2 -p 12 --very-sensitive --rg-id K562_RNA-seq_10 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/tmpui_j0o7u -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam` (67729,67735,67736)
<pre>
62592060 reads; of these:
  62592060 (100.00%) were unpaired; of these:
    1837492 (2.94%) aligned 0 times
    44995618 (71.89%) aligned exactly 1 time
    15758950 (25.18%) aligned >1 times
97.06% overall alignment rate
[bam_sort_core] merging from 20 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:45:05. Running peak memory: 4.589GB.  
  PID: 67729;	Command: bowtie2;	Return code: 0;	Memory used: 3.72GB  
  PID: 67735;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 67736;	Command: samtools;	Return code: 0;	Memory used: 0.889GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam` (91421)
<pre>
</pre>
Command completed. Elapsed time: 0:02:10. Running peak memory: 4.589GB.  
  PID: 91421;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	60754568	PEPPRO	_RES_

> `QC_filtered_reads`	7195486	PEPPRO	_RES_

> `Aligned_reads`	53559082	PEPPRO	_RES_

> `Alignment_rate`	78.28	PEPPRO	_RES_

> `Total_efficiency`	76.51	PEPPRO	_RES_

> `Read_depth`	5.9	PEPPRO	_RES_

### Compress all unmapped read files (06-11 18:41:44) elapsed: 4202.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq` (94676)
<pre>
</pre>
Command completed. Elapsed time: 0:02:17. Running peak memory: 4.589GB.  
  PID: 94676;	Command: pigz;	Return code: 0;	Memory used: 0.011GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam` (94860)
<pre>
</pre>
Command completed. Elapsed time: 0:01:04. Running peak memory: 4.589GB.  
  PID: 94860;	Command: samtools;	Return code: 0;	Memory used: 0.017GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	1579435	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam` (95156)
<pre>
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 4.589GB.  
  PID: 95156;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/chr_sizes.bed` (95253,95254,95255,95256)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.589GB.  
  PID: 95255;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 95253;	Command: samtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 95256;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 95254;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_noMT.bam` (95258)
<pre>
</pre>
Command completed. Elapsed time: 0:01:14. Running peak memory: 4.589GB.  
  PID: 95258;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam` (95367)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.589GB.  
  PID: 95367;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam` (95368)
<pre>
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 4.589GB.  
  PID: 95368;	Command: samtools;	Return code: 0;	Memory used: 0.015GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-11 18:50:58) elapsed: 555.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv` (95810)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/tmp_K562_RNA-seq_10_sort_sqa2liyh'
Processing with 12 cores...
Discarding 97 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270725v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 98 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:01:24. Running peak memory: 4.589GB.  
  PID: 95810;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.428GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv`

> `NRF`	0.77	PEPPRO	_RES_

> `PBC1`	0.89	PEPPRO	_RES_

> `PBC2`	11.77	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_unmap.bam`  

> `samtools view -b -@ 12 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_unmap.bam` (95913)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 4.589GB.  
  PID: 95913;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_temp.bam`

> `Unmapped_reads`	1837492	PEPPRO	_RES_

### Split BAM by strand (06-11 18:52:45) elapsed: 107.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam` (95962)
<pre>
</pre>
Got SIGTERM. Failing gracefully... (06-11 18:55:14) elapsed: 149.0 _TIME_
Child process 95962 (samtools) terminated after 1 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/recover.lock.aligned_hg38__K562_RNA-seq_10_plus.bam
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/recover.lock.aligned_hg38__K562_RNA-seq_10_minus.bam

### Pipeline failed at:  (06-11 18:55:15) elapsed: 1.0 _TIME_

Total time: 1:46:53
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_10 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_10pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba26-34c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/
*  Pipeline started at:   (06-11 19:08:10) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_10pctRNA.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_10`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_10pctRNA.fastq.gz

> `File_mb`	5660.06	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 19:08:11) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz']

### FASTQ processing:  (06-11 19:08:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-11 19:08:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-11 19:08:11) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 19:08:11) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq

### Map to genome (06-11 19:08:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam`  

### Compress all unmapped read files (06-11 19:08:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-11 19:08:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_unmap.bam`  

### Split BAM by strand (06-11 19:08:11) elapsed: 0.0 _TIME_

Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/lock.aligned_hg38__K562_RNA-seq_10_plus.bam
Overwriting target...
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/lock.aligned_hg38__K562_RNA-seq_10_minus.bam
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam` (199879)
<pre>
</pre>
Command completed. Elapsed time: 0:04:18. Running peak memory: 0.006GB.  
  PID: 199879;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam` (200545)
<pre>
</pre>
Command completed. Elapsed time: 0:04:07. Running peak memory: 0.006GB.  
  PID: 200545;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-11 19:16:36) elapsed: 505.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (201397)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.006GB.  
  PID: 201397;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_plus_TssEnrichment.txt` (201398)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 0.535GB.  
  PID: 201398;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.535GB


> `TSS_coding_score`	15.7	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_minus_TssEnrichment.txt` (201436)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 0.546GB.  
  PID: 201436;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.546GB


> `TSS_non-coding_score`	5.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_minus_TssEnrichment.txt` (201470)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.546GB.  
  PID: 201470;	Command: Rscript;	Return code: 0;	Memory used: 0.285GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_10_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_10_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt` (201494,201495,201496,201497)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.546GB.  
  PID: 201494;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 201496;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 201495;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 201497;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt` (201499)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.546GB.  
  PID: 201499;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-11 19:17:02) elapsed: 26.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_tss.bed` (201501,201502)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.546GB.  
  PID: 201501;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 201502;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_gene_body.bed` (201505,201506)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.546GB.  
  PID: 201505;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 201506;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSS_density.bed` (201509,201510,201511,201512)
<pre>
</pre>
Command completed. Elapsed time: 0:01:20. Running peak memory: 0.546GB.  
  PID: 201510;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 201512;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 201509;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB  
  PID: 201511;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_gene_body_density.bed` (201698,201699,201700)
<pre>
</pre>
Command completed. Elapsed time: 0:02:03. Running peak memory: 0.546GB.  
  PID: 201699;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 201698;	Command: bedtools;	Return code: 0;	Memory used: 0.155GB  
  PID: 201700;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/tmphkr_pqiv` (202108,202109,202110)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.546GB.  
  PID: 202108;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 202110;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 202109;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/tmphkr_pqiv | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-11 19:20:28) elapsed: 206.0 _TIME_

Total time: 0:12:19
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/tmphkr_pqiv | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_10 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_10pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aw29-23b
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/
*  Pipeline started at:   (06-14 21:11:17) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_10pctRNA.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_10`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_10pctRNA.fastq.gz

> `File_mb`	5660.06	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-14 21:11:17) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/K562_RNA-seq_10_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/K562_RNA-seq_10.fastq.gz']

### FASTQ processing:  (06-14 21:11:17) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-14 21:11:17) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_10_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-14 21:11:17) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:11:17) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq

### Map to genome (06-14 21:11:17) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam`  

### Compress all unmapped read files (06-14 21:11:17) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/prealignments/K562_RNA-seq_10_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:11:17) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_unmap.bam`  

### Split BAM by strand (06-14 21:11:17) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam`  

### Calculate TSS enrichment (06-14 21:11:17) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/K562_RNA-seq_10_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_10_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:11:17) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/tmp1o6x5gfq` (425686,425687,425688)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 425686;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 425688;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 425687;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/tmp1o6x5gfq | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0420842) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/tmp1o6x5gfq > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed` (425694)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 425694;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	9.88	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed` (425699)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.22GB.  
  PID: 425699;	Command: Rscript;	Return code: 0;	Memory used: 0.22GB

> `Pause index`	QC_hg38/K562_RNA-seq_10_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_10_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_pause_index.bed` (425722)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.22GB.  
  PID: 425722;	Command: pigz;	Return code: 0;	Memory used: 0.001GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:11:24) elapsed: 7.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam`
53559082 19292029

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam`
53559082 18473092

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_gene_sort.bed` (429552,429570)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.22GB.  
  PID: 429552;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 429570;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_gene_coverage.bed` (429911)
<pre>
</pre>
Command completed. Elapsed time: 0:02:10. Running peak memory: 0.22GB.  
  PID: 429911;	Command: bedtools;	Return code: 0;	Memory used: 0.154GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed.gz` (6234)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.22GB.  
  PID: 6234;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed` (6235)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.22GB.  
  PID: 6235;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:15:15) elapsed: 231.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/raw/hg38_annotations.bed` (6244)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.22GB.  
  PID: 6244;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer_sort.bed` (6249,6250,6251,6252)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.22GB.  
  PID: 6249;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 6250;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 6252;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB  
  PID: 6251;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_plus_coverage.bed` (6255)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 0.22GB.  
  PID: 6255;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_minus_coverage.bed` (6289)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 0.22GB.  
  PID: 6289;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_sort.bed` (6327,6328,6329,6330)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.22GB.  
  PID: 6327;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 6328;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 6330;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 6329;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_plus_coverage.bed` (6332)
<pre>
</pre>
Command completed. Elapsed time: 0:00:48. Running peak memory: 0.22GB.  
  PID: 6332;	Command: bedtools;	Return code: 0;	Memory used: 0.091GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_minus_coverage.bed` (23860)
<pre>
</pre>
Command completed. Elapsed time: 0:00:47. Running peak memory: 0.22GB.  
  PID: 23860;	Command: bedtools;	Return code: 0;	Memory used: 0.07GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region"` (39929)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.22GB.  
  PID: 39929;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed` (40211,40212,40213,40214)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.22GB.  
  PID: 40211;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 40213;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 40212;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 40214;	Command: bedtools;	Return code: 0;	Memory used: 0.05GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_plus_coverage.bed` (40498)
<pre>
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 0.22GB.  
  PID: 40498;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_minus_coverage.bed` (45163)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 0.22GB.  
  PID: 45163;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR"` (45248)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.22GB.  
  PID: 45248;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR_sort.bed` (45250,45251,45252,45253)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.22GB.  
  PID: 45250;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 45251;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 45253;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 45252;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_plus_coverage.bed` (45262)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 0.22GB.  
  PID: 45262;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_minus_coverage.bed` (45601)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 0.22GB.  
  PID: 45601;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR"` (45636)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.22GB.  
  PID: 45636;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR_sort.bed` (45637,45638,45639,45640)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.22GB.  
  PID: 45637;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 45638;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 45640;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 45639;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_plus_coverage.bed` (45643)
<pre>
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 0.22GB.  
  PID: 45643;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_minus_coverage.bed` (45682)
<pre>
</pre>
Command completed. Elapsed time: 0:00:43. Running peak memory: 0.22GB.  
  PID: 45682;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon_sort.bed` (51004,51005,51006,51010)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.22GB.  
  PID: 51004;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 51005;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 51010;	Command: bedtools;	Return code: 0;	Memory used: 0.173GB  
  PID: 51006;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_plus_coverage.bed` (53399)
<pre>
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 0.22GB.  
  PID: 53399;	Command: bedtools;	Return code: 0;	Memory used: 0.081GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_minus_coverage.bed` (78822)
<pre>
</pre>
Command completed. Elapsed time: 0:00:45. Running peak memory: 0.22GB.  
  PID: 78822;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron_sort.bed` (84365,84366,84367,84368)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.22GB.  
  PID: 84365;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 84367;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 84366;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 84368;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_plus_coverage.bed` (84372)
<pre>
</pre>
Command completed. Elapsed time: 0:00:56. Running peak memory: 0.22GB.  
  PID: 84372;	Command: bedtools;	Return code: 0;	Memory used: 0.069GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_minus_coverage.bed` (86118)
<pre>
</pre>
Command completed. Elapsed time: 0:00:53. Running peak memory: 0.22GB.  
  PID: 86118;	Command: bedtools;	Return code: 0;	Memory used: 0.069GB


### Plot cFRiF/FRiF (06-14 21:25:57) elapsed: 643.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_RNA-seq_10 -z 3099922541 -n 26481130 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_plus_coverage.bed` (107631)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:44. Running peak memory: 0.436GB.  
  PID: 107631;	Command: Rscript;	Return code: 0;	Memory used: 0.436GB

> `cFRiF`	QC_hg38/K562_RNA-seq_10_cFRiF.pdf	cFRiF	QC_hg38/K562_RNA-seq_10_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_RNA-seq_10 -z 3099922541 -n 26481130 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_Intron_plus_coverage.bed` (123186)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 0.436GB.  
  PID: 123186;	Command: Rscript;	Return code: 0;	Memory used: 0.436GB

> `FRiF`	QC_hg38/K562_RNA-seq_10_FRiF.pdf	FRiF	QC_hg38/K562_RNA-seq_10_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 21:27:16) elapsed: 79.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_exons_sort.bed` (123227,123228)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.436GB.  
  PID: 123228;	Command: bedtools;	Return code: 0;	Memory used: 0.095GB  
  PID: 123227;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_introns_sort.bed` (123238,123239,123240)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.436GB.  
  PID: 123238;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 123240;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 123239;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_coverage.bed` (123258)
<pre>
</pre>
Command completed. Elapsed time: 0:01:25. Running peak memory: 0.436GB.  
  PID: 123258;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_coverage.bed` (123361)
<pre>
</pre>
Command completed. Elapsed time: 0:01:46. Running peak memory: 0.436GB.  
  PID: 123361;	Command: bedtools;	Return code: 0;	Memory used: 0.09GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/53.559082)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_rpkm.bed` (142703,142704,142705)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.436GB.  
  PID: 142703;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 142705;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 142704;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/53.559082)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_rpkm.bed` (142857,142861,142865)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.436GB.  
  PID: 142857;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 142865;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 142861;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed` (143040,143044,143064)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.436GB.  
  PID: 143040;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 143064;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 143044;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.97	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed --annotate` (143151)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.436GB.  
  PID: 143151;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_10_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_10_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/QC_hg38/K562_RNA-seq_10_exon_intron_ratios.bed` (144522)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.436GB.  
  PID: 144522;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-14 21:30:44) elapsed: 208.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam` (144545)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 0.436GB.  
  PID: 144545;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 53559082.0` (147483)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_10_plus_cuttrace_1e_y28uq'
Processing with 4 cores...
stdin is empty of data
Discarding 107 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270732v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 88 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 88 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_exact_body_0-mer.bw'
Merging 88 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:58. Running peak memory: 3.824GB.  
  PID: 147483;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.824GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam` (164216)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.824GB.  
  PID: 164216;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 53559082.0` (164349)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/aligned_hg38/K562_RNA-seq_10_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_10_minus_cuttrace_hl_k5fsm'
Processing with 4 cores...
stdin is empty of data
stdin is empty of data
Discarding 109 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1']
Keeping 86 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_exact_body_0-mer.bw'
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_10/signal_hg38/K562_RNA-seq_10_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:27. Running peak memory: 3.824GB.  
  PID: 164349;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.645GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:37:46
*  Total elapsed time (all runs):  3:01:29
*         Peak memory (this run):  3.8243 GB
*        Pipeline completed time: 2020-06-14 21:49:02
