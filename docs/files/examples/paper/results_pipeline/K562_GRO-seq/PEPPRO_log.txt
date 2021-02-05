### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_GRO-seq --genome hg38 --input /project/shefflab/data//sra_fastq/SRR1552484.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol GRO --umi-len 0 --scale --prealignments human_rDNA --dirty`
*         Compute host:  udc-aj38-17c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/
*  Pipeline started at:   (06-11 17:29:51) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//sra_fastq/SRR1552484.fastq.gz']`
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
*           `protocol`:  `GRO`
*            `recover`:  `False`
*        `sample_name`:  `K562_GRO-seq`
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

Local input file: /project/shefflab/data//sra_fastq/SRR1552484.fastq.gz

> `File_mb`	1322.25	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected GRO input

### Merge/link and fastq conversion:  (06-11 17:29:53) elapsed: 2.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz`  

> `ln -sf /project/shefflab/data//sra_fastq/SRR1552484.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz` (347126)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 347126;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1.fastq` (347127)
<pre>
</pre>
Command completed. Elapsed time: 0:01:06. Running peak memory: 0.002GB.  
  PID: 347127;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	30363736	PEPPRO	_RES_

> `Fastq_reads`	30363736	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz']

### FASTQ processing:  (06-11 17:31:22) elapsed: 89.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1_processed.fastq`  

> `(cutadapt -j 12 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt` (347505)
<pre>
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 2.13GB.  
  PID: 347505;	Command: cutadapt;	Return code: 0;	Memory used: 2.13GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1_noadap.fastq | seqtk seq -L 2 - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1_processed.fastq` (347576,347578)
<pre>
</pre>
Command completed. Elapsed time: 0:01:00. Running peak memory: 2.13GB.  
  PID: 347576;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 347578;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	22296806.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	1041369.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	3.4296	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastqc/K562_GRO-seq_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (347706)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.13GB.  
  PID: 347706;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	29322367	PEPPRO	_RES_

> `Trim_loss_rate`	3.43	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1_processed.fastq` (347737)
<pre>
Started analysis of K562_GRO-seq_R1_processed.fastq
Approx 5% complete for K562_GRO-seq_R1_processed.fastq
Approx 10% complete for K562_GRO-seq_R1_processed.fastq
Approx 15% complete for K562_GRO-seq_R1_processed.fastq
Approx 20% complete for K562_GRO-seq_R1_processed.fastq
Approx 25% complete for K562_GRO-seq_R1_processed.fastq
Approx 30% complete for K562_GRO-seq_R1_processed.fastq
Approx 35% complete for K562_GRO-seq_R1_processed.fastq
Approx 40% complete for K562_GRO-seq_R1_processed.fastq
Approx 45% complete for K562_GRO-seq_R1_processed.fastq
Approx 50% complete for K562_GRO-seq_R1_processed.fastq
Approx 55% complete for K562_GRO-seq_R1_processed.fastq
Approx 60% complete for K562_GRO-seq_R1_processed.fastq
Approx 65% complete for K562_GRO-seq_R1_processed.fastq
Approx 70% complete for K562_GRO-seq_R1_processed.fastq
Approx 75% complete for K562_GRO-seq_R1_processed.fastq
Approx 80% complete for K562_GRO-seq_R1_processed.fastq
Approx 85% complete for K562_GRO-seq_R1_processed.fastq
Approx 90% complete for K562_GRO-seq_R1_processed.fastq
Approx 95% complete for K562_GRO-seq_R1_processed.fastq
Analysis complete for K562_GRO-seq_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:29. Running peak memory: 2.13GB.  
  PID: 347737;	Command: fastqc;	Return code: 0;	Memory used: 0.165GB

> `FastQC report r1`	fastqc/K562_GRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/processed_R1.flag` (348047)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.13GB.  
  PID: 348047;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Plot adapter insertion distribution (06-11 17:35:38) elapsed: 256.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt` (348048)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 2.13GB.  
  PID: 348048;	Command: Rscript;	Return code: 0;	Memory used: 0.112GB

> `Adapter insertion distribution`	cutadapt/K562_GRO-seq_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_GRO-seq_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	22	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 17:35:49) elapsed: 11.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.7189	PEPPRO	_RES_

### Prealignments (06-11 17:35:49) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 17:35:49) elapsed: 0.0 _TIME_


> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_GRO-seq -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/prealignments/K562_GRO-seq_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
29322367 reads; of these:
  29322367 (100.00%) were unpaired; of these:
    25012896 (85.30%) aligned 0 times
    4309471 (14.70%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
14.70% overall alignment rate

> `Aligned_reads_human_rDNA`	4309471.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	14.7	PEPPRO	_RES_

### Map to genome (06-11 17:39:22) elapsed: 213.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam`  

> `bowtie2 -p 12 --very-sensitive --rg-id K562_GRO-seq -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/prealignments/K562_GRO-seq_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/tmpuee4zdhj -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_temp.bam` (348363,348364,348365)
<pre>
25012896 reads; of these:
  25012896 (100.00%) were unpaired; of these:
    2113226 (8.45%) aligned 0 times
    15215880 (60.83%) aligned exactly 1 time
    7683790 (30.72%) aligned >1 times
91.55% overall alignment rate
[bam_sort_core] merging from 6 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:10:14. Running peak memory: 3.668GB.  
  PID: 348363;	Command: bowtie2;	Return code: 0;	Memory used: 3.668GB  
  PID: 348364;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 348365;	Command: samtools;	Return code: 0;	Memory used: 0.902GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam` (349569)
<pre>
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 3.668GB.  
  PID: 349569;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	22899670	PEPPRO	_RES_

> `QC_filtered_reads`	4468262	PEPPRO	_RES_

> `Aligned_reads`	18431408	PEPPRO	_RES_

> `Alignment_rate`	62.86	PEPPRO	_RES_

> `Total_efficiency`	60.7	PEPPRO	_RES_

> `Read_depth`	2.65	PEPPRO	_RES_

### Compress all unmapped read files (06-11 18:00:14) elapsed: 1252.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/prealignments/K562_GRO-seq_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/prealignments/K562_GRO-seq_human_rDNA_unmap.fq` (351183)
<pre>
</pre>
Command completed. Elapsed time: 0:00:47. Running peak memory: 3.668GB.  
  PID: 351183;	Command: pigz;	Return code: 0;	Memory used: 0.013GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_temp.bam` (351262)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 3.668GB.  
  PID: 351262;	Command: samtools;	Return code: 0;	Memory used: 0.017GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	843399	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam` (351285)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.668GB.  
  PID: 351285;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/chr_sizes.bed` (351301,351302,351303,351304)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.668GB.  
  PID: 351301;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 351303;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 351302;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 351304;	Command: grep;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_noMT.bam` (351307)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 3.668GB.  
  PID: 351307;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam` (351355)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.668GB.  
  PID: 351355;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam` (351356)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.668GB.  
  PID: 351356;	Command: samtools;	Return code: 0;	Memory used: 0.014GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	50	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-11 18:03:14) elapsed: 180.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_bamQC.tsv` (351434)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/tmp_K562_GRO-seq_sort_gbhq8men'
Processing with 12 cores...
Discarding 91 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1']
Keeping 104 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270511v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270374v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.668GB.  
  PID: 351434;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.057GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_bamQC.tsv`

> `NRF`	0.7	PEPPRO	_RES_

> `PBC1`	0.86	PEPPRO	_RES_

> `PBC2`	8.69	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_unmap.bam`  

> `samtools view -b -@ 12 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_unmap.bam` (351500)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.668GB.  
  PID: 351500;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_temp.bam`

> `Unmapped_reads`	2113226	PEPPRO	_RES_

### Split BAM by strand (06-11 18:03:58) elapsed: 44.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam` (351538)
<pre>
</pre>
Command completed. Elapsed time: 0:01:06. Running peak memory: 3.668GB.  
  PID: 351538;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam` (351800)
<pre>
</pre>
Command completed. Elapsed time: 0:01:03. Running peak memory: 3.668GB.  
  PID: 351800;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-11 18:06:07) elapsed: 129.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (351866)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.668GB.  
  PID: 351866;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_plus_TssEnrichment.txt` (351867)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.668GB.  
  PID: 351867;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.599GB


> `TSS_coding_score`	23.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_minus_TssEnrichment.txt` (351903)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.668GB.  
  PID: 351903;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.746GB


> `TSS_non-coding_score`	8.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_minus_TssEnrichment.txt` (351934)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.668GB.  
  PID: 351934;	Command: Rscript;	Return code: 0;	Memory used: 0.102GB

> `TSS enrichment`	QC_hg38/K562_GRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_GRO-seq_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt` (351961,351962,351963,351964)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.668GB.  
  PID: 351961;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 351963;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 351962;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 351964;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_keep.txt` (351967)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.668GB.  
  PID: 351967;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-11 18:06:32) elapsed: 25.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_tss.bed` (351969,351970)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.668GB.  
  PID: 351969;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 351970;	Command: bedtools;	Return code: 0;	Memory used: 0.094GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed` (351974,351975)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.668GB.  
  PID: 351974;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 351975;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_TSS_density.bed` (351979,351980,351981,351982)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.668GB.  
  PID: 351979;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB  
  PID: 351981;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 351980;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 351982;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_gene_body_density.bed` (352009,352010,352012)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 3.668GB.  
  PID: 352012;	Command: sort;	Return code: 0;	Memory used: 0.007GB  
  PID: 352009;	Command: bedtools;	Return code: 0;	Memory used: 0.051GB  
  PID: 352010;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed`  
awk: cmd. line:1: { if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}
awk: cmd. line:1:             ^ syntax error
awk: cmd. line:1: { if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}
awk: cmd. line:1:                                                                                                                      ^ syntax error

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/tmpxwvopkcb` (352054,352055,352056)
<pre>
</pre>
Got SIGTERM. Failing gracefully... (06-11 18:55:15) elapsed: 2923.0 _TIME_
Child process 352054 (join) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/recover.lock.QC_hg38__K562_GRO-seq_pause_index.bed

### Pipeline failed at:  (06-11 18:55:15) elapsed: 0.0 _TIME_

Total time: 1:25:24
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_GRO-seq --genome hg38 --input /project/shefflab/data//sra_fastq/SRR1552484.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol GRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba26-16
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/
*  Pipeline started at:   (06-11 19:10:15) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//sra_fastq/SRR1552484.fastq.gz']`
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
*           `protocol`:  `GRO`
*            `recover`:  `True`
*        `sample_name`:  `K562_GRO-seq`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//sra_fastq/SRR1552484.fastq.gz

> `File_mb`	1322.25	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected GRO input

### Merge/link and fastq conversion:  (06-11 19:10:15) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz']

### FASTQ processing:  (06-11 19:10:15) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-11 19:10:15) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_GRO-seq_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_GRO-seq_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-11 19:10:15) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 19:10:15) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/prealignments/K562_GRO-seq_human_rDNA_unmap.fq

### Map to genome (06-11 19:10:15) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam`  

### Compress all unmapped read files (06-11 19:10:15) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/prealignments/K562_GRO-seq_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-11 19:10:15) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_unmap.bam`  

### Split BAM by strand (06-11 19:10:15) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam`  

### Calculate TSS enrichment (06-11 19:10:15) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/K562_GRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_GRO-seq_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-11 19:10:15) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_gene_body_density.bed`  
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/lock.QC_hg38__K562_GRO-seq_pause_index.bed
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/tmpnu6wpr_c` (71328,71329,71330)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.006GB.  
  PID: 71328;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 71330;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 71329;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/tmpnu6wpr_c | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-11 19:10:16) elapsed: 0.0 _TIME_

Total time: 0:00:01
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/tmpnu6wpr_c | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_GRO-seq --genome hg38 --input /project/shefflab/data//sra_fastq/SRR1552484.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol GRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aj37-18c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/
*  Pipeline started at:   (06-14 21:14:10) elapsed: 9.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//sra_fastq/SRR1552484.fastq.gz']`
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
*           `protocol`:  `GRO`
*            `recover`:  `True`
*        `sample_name`:  `K562_GRO-seq`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//sra_fastq/SRR1552484.fastq.gz

> `File_mb`	1322.25	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected GRO input

### Merge/link and fastq conversion:  (06-14 21:14:11) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz']

### FASTQ processing:  (06-14 21:14:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-14 21:14:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_GRO-seq_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_GRO-seq_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-14 21:14:11) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:14:11) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/prealignments/K562_GRO-seq_human_rDNA_unmap.fq

### Map to genome (06-14 21:14:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam`  

### Compress all unmapped read files (06-14 21:14:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/prealignments/K562_GRO-seq_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:14:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_unmap.bam`  

### Split BAM by strand (06-14 21:14:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam`  

### Calculate TSS enrichment (06-14 21:14:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/K562_GRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_GRO-seq_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:14:11) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/tmp4yy_drxc` (428183,428191,428201)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 428183;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 428201;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 428191;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/tmp4yy_drxc | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0172824) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/tmp4yy_drxc > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed` (428280)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 428280;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.74	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed` (428320)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.287GB.  
  PID: 428320;	Command: Rscript;	Return code: 0;	Memory used: 0.287GB

> `Pause index`	QC_hg38/K562_GRO-seq_pause_index.pdf	Pause index	QC_hg38/K562_GRO-seq_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed` (429780)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.287GB.  
  PID: 429780;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:14:18) elapsed: 6.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam`
18431408 6290701

> `Plus_FRiP`	0.34	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam`
18431408 5983864

> `Minus_FRiP`	0.32	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_gene_sort.bed` (430842,430843)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.287GB.  
  PID: 430842;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 430843;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_gene_coverage.bed` (430845)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 0.287GB.  
  PID: 430845;	Command: bedtools;	Return code: 0;	Memory used: 0.074GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/hg38_annotations.bed.gz` (431095)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.287GB.  
  PID: 431095;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/hg38_annotations.bed` (431096)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.287GB.  
  PID: 431096;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:15:17) elapsed: 60.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/raw/hg38_annotations.bed` (431105)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.287GB.  
  PID: 431105;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Enhancer_sort.bed` (431108,431109,431111,431113)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.287GB.  
  PID: 431108;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 431109;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 431113;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB  
  PID: 431111;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Enhancer_plus_coverage.bed` (431115)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 0.287GB.  
  PID: 431115;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Enhancer_minus_coverage.bed` (431133)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 0.287GB.  
  PID: 431133;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_sort.bed` (431145,431146,431147,431148)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.287GB.  
  PID: 431145;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 431146;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 431148;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 431147;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_plus_coverage.bed` (431150)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 0.287GB.  
  PID: 431150;	Command: bedtools;	Return code: 0;	Memory used: 0.025GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_minus_coverage.bed` (431171)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 0.287GB.  
  PID: 431171;	Command: bedtools;	Return code: 0;	Memory used: 0.027GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_Flanking_Region"` (431197)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.287GB.  
  PID: 431197;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed` (431198,431199,431200,431201)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.287GB.  
  PID: 431198;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 431200;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 431199;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 431201;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (431204)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 0.287GB.  
  PID: 431204;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_Flanking_Region_minus_coverage.bed` (431228)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 0.287GB.  
  PID: 431228;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/5_UTR"` (431248)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.287GB.  
  PID: 431248;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/5_UTR_sort.bed` (431249,431250,431251,431252)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.287GB.  
  PID: 431249;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 431250;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 431252;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 431251;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_5_UTR_plus_coverage.bed` (431255)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 0.287GB.  
  PID: 431255;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_5_UTR_minus_coverage.bed` (431268)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 0.287GB.  
  PID: 431268;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/3_UTR"` (431285)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.287GB.  
  PID: 431285;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/3_UTR_sort.bed` (431286,431287,431288,431289)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.287GB.  
  PID: 431286;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 431287;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 431289;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 431288;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_3_UTR_plus_coverage.bed` (431291)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 0.287GB.  
  PID: 431291;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_3_UTR_minus_coverage.bed` (431349)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 0.287GB.  
  PID: 431349;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Exon_sort.bed` (436564,436572,436585,436590)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.287GB.  
  PID: 436564;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 436572;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 436590;	Command: bedtools;	Return code: 0;	Memory used: 0.172GB  
  PID: 436585;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Exon_plus_coverage.bed` (439659)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 0.287GB.  
  PID: 439659;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Exon_minus_coverage.bed` (456903)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 0.287GB.  
  PID: 456903;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Intron_sort.bed` (14009,14016,14019,14023)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.287GB.  
  PID: 14009;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 14019;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 14016;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 14023;	Command: bedtools;	Return code: 0;	Memory used: 0.078GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Intron_plus_coverage.bed` (15300)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 0.287GB.  
  PID: 15300;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Intron_minus_coverage.bed` (17150)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 0.287GB.  
  PID: 17150;	Command: bedtools;	Return code: 0;	Memory used: 0.028GB


### Plot cFRiF/FRiF (06-14 21:18:00) elapsed: 163.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_GRO-seq -z 3099922541 -n 9058650 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Intron_plus_coverage.bed` (21478)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 0.462GB.  
  PID: 21478;	Command: Rscript;	Return code: 0;	Memory used: 0.462GB

> `cFRiF`	QC_hg38/K562_GRO-seq_cFRiF.pdf	cFRiF	QC_hg38/K562_GRO-seq_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_GRO-seq -z 3099922541 -n 9058650 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Intron_plus_coverage.bed` (53701)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 0.507GB.  
  PID: 53701;	Command: Rscript;	Return code: 0;	Memory used: 0.507GB

> `FRiF`	QC_hg38/K562_GRO-seq_FRiF.pdf	FRiF	QC_hg38/K562_GRO-seq_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 21:19:00) elapsed: 61.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_exons_sort.bed` (61345,61346)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.507GB.  
  PID: 61345;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 61346;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_introns_sort.bed` (61352,61353,61354)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.507GB.  
  PID: 61352;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 61354;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 61353;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exons_coverage.bed` (61366)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 0.507GB.  
  PID: 61366;	Command: bedtools;	Return code: 0;	Memory used: 0.047GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_introns_coverage.bed` (61487)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 0.507GB.  
  PID: 61487;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.431408)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exons_rpkm.bed` (61511,61512,61513)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.507GB.  
  PID: 61511;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 61513;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 61512;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.431408)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_introns_rpkm.bed` (61516,61517,61518)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.507GB.  
  PID: 61516;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 61518;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 61517;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exon_intron_ratios.bed` (61521,61522,61523)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.507GB.  
  PID: 61521;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 61523;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 61522;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.33	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exon_intron_ratios.bed --annotate` (61529)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.507GB.  
  PID: 61529;	Command: Rscript;	Return code: 0;	Memory used: 0.323GB

> `mRNA contamination`	QC_hg38/K562_GRO-seq_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_GRO-seq_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exon_intron_ratios.bed` (61559)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.507GB.  
  PID: 61559;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (06-14 21:19:59) elapsed: 58.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam` (61567)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.507GB.  
  PID: 61567;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_plus_smooth_body_0-mer.bw -p 8 --variable-step --scale 18431408.0` (61790)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam'
Temporary files will be stored in: 'tmp_K562_GRO-seq_plus_cuttrace_pp7ztq7s'
Processing with 4 cores...
stdin is empty of data
stdin is empty of data
Discarding 100 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr14_KI270723v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1']
Keeping 95 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270511v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270374v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 95 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_plus_exact_body_0-mer.bw'
Merging 95 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:24. Running peak memory: 2.62GB.  
  PID: 61790;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.62GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam` (239914)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 2.62GB.  
  PID: 239914;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_minus_smooth_body_0-mer.bw -p 8 --variable-step --scale 18431408.0` (239935)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam'
Temporary files will be stored in: 'tmp_K562_GRO-seq_minus_cuttrace_1axk6ebz'
Processing with 4 cores...
stdin is empty of data
Discarding 111 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr11_KI270721v1_random', 'chr17_KI270730v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 84 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_minus_exact_body_0-mer.bw'
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:16. Running peak memory: 2.732GB.  
  PID: 239935;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.732GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:20:51
*  Total elapsed time (all runs):  0:58:31
*         Peak memory (this run):  2.7315 GB
*        Pipeline completed time: 2020-06-14 21:34:52
