### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq --genome hg38 --input /project/shefflab/data//sra_fastq/SRR1554311.fastq.gz /project/shefflab/data//sra_fastq/SRR1554312.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/results_pipeline -P 24 -M 24000 --protocol PRO --prealignments human_rDNA`
*         Compute host:  udc-aj37-15c0
*          Working dir:  /sfs/qumulo/qproject/shefflab/processed/peppro/paper/submission
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/
*  Pipeline started at:   (05-27 08:15:11) elapsed: 12.0 _TIME_

### Version log:

*       Python version:  3.6.6
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.9.7
*        Pipeline hash:  70861b332c5d9a7ad285ab85db3c0f2a94a3f099
*      Pipeline branch:  * master
*        Pipeline date:  2020-05-26 16:17:03 -0400

### Arguments passed to pipeline:

*           `TSS_name`:  `None`
*            `adapter`:  `cutadapt`
*          `anno_name`:  `None`
*         `complexity`:  `False`
*        `config_file`:  `peppro.yaml`
*              `cores`:  `24`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `False`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//sra_fastq/SRR1554311.fastq.gz', '/project/shefflab/data//sra_fastq/SRR1554312.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `24000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `False`
*        `sample_name`:  `K562_PRO-seq`
*              `scale`:  `False`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `SINGLE`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `0`
*          `verbosity`:  `None`

----------------------------------------

Local input file: /project/shefflab/data//sra_fastq/SRR1554311.fastq.gz

> `File_mb`	38578.5	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (05-27 08:15:11) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz`  

> `cat /project/shefflab/data//sra_fastq/SRR1554311.fastq.gz /project/shefflab/data//sra_fastq/SRR1554312.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz` (70931)
<pre>
</pre>
Command completed. Elapsed time: 0:03:41. Running peak memory: 0.003GB.  
  PID: 70931;	Command: cat;	Return code: 0;	Memory used: 0.003GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1.fastq`  

> `pigz -f -p 24 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1.fastq` (71154)
<pre>
</pre>
Command completed. Elapsed time: 0:18:00. Running peak memory: 0.003GB.  
  PID: 71154;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	496581677	PEPPRO	_RES_

> `Fastq_reads`	496581677	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/K562_PRO-seq.merged.fastq.gz']

### FASTQ processing:  (05-27 08:50:05) elapsed: 2094.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_processed.fastq`  

> `(cutadapt -j 24 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt` (84332)
<pre>
</pre>
Command completed. Elapsed time: 0:13:47. Running peak memory: 0.748GB.  
  PID: 84332;	Command: cutadapt;	Return code: 0;	Memory used: 0.748GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_processed.fastq` (87891,87892)
<pre>
</pre>
Command completed. Elapsed time: 0:09:51. Running peak memory: 0.748GB.  
  PID: 87891;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 87892;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	423059183.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Uninformative_adapter_reads`	0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	0.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastqc/K562_PRO-seq_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (89531)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.748GB.  
  PID: 89531;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	484143694	PEPPRO	_RES_

> `Trim_loss_rate`	2.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_processed.fastq` (89927)
<pre>
Started analysis of K562_PRO-seq_R1_processed.fastq
Approx 5% complete for K562_PRO-seq_R1_processed.fastq
Approx 10% complete for K562_PRO-seq_R1_processed.fastq
Approx 15% complete for K562_PRO-seq_R1_processed.fastq
Approx 20% complete for K562_PRO-seq_R1_processed.fastq
Approx 25% complete for K562_PRO-seq_R1_processed.fastq
Approx 30% complete for K562_PRO-seq_R1_processed.fastq
Approx 35% complete for K562_PRO-seq_R1_processed.fastq
Approx 40% complete for K562_PRO-seq_R1_processed.fastq
Approx 45% complete for K562_PRO-seq_R1_processed.fastq
Approx 50% complete for K562_PRO-seq_R1_processed.fastq
Approx 55% complete for K562_PRO-seq_R1_processed.fastq
Approx 60% complete for K562_PRO-seq_R1_processed.fastq
Approx 65% complete for K562_PRO-seq_R1_processed.fastq
Approx 70% complete for K562_PRO-seq_R1_processed.fastq
Approx 75% complete for K562_PRO-seq_R1_processed.fastq
Approx 80% complete for K562_PRO-seq_R1_processed.fastq
Approx 85% complete for K562_PRO-seq_R1_processed.fastq
Approx 90% complete for K562_PRO-seq_R1_processed.fastq
Approx 95% complete for K562_PRO-seq_R1_processed.fastq
Analysis complete for K562_PRO-seq_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:30:25. Running peak memory: 0.748GB.  
  PID: 89927;	Command: fastqc;	Return code: 0;	Memory used: 0.306GB

> `FastQC report r1`	fastqc/K562_PRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/processed_R1.flag` (93279)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.748GB.  
  PID: 93279;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (05-27 09:51:41) elapsed: 3696.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt` (93285)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.748GB.  
  PID: 93285;	Command: Rscript;	Return code: 0;	Memory used: 0.201GB

> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (05-27 09:51:48) elapsed: 7.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/cutadapt/K562_PRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2312	PEPPRO	_RES_

### Prealignments (05-27 09:51:48) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (05-27 09:51:48) elapsed: 0.0 _TIME_


> `(bowtie2 -p 24 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_PRO-seq -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/fastq/K562_PRO-seq_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/prealignments/K562_PRO-seq_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
484143694 reads; of these:
  484143694 (100.00%) were unpaired; of these:
    439463330 (90.77%) aligned 0 times
    44680364 (9.23%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.23% overall alignment rate

> `Aligned_reads_human_rDNA`	44680364.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.23	PEPPRO	_RES_

### Map to genome (05-27 10:49:38) elapsed: 3471.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam`  

> `bowtie2 -p 24 --very-sensitive --rg-id K562_PRO-seq -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/prealignments/K562_PRO-seq_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/tmpvvymivrz -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam` (99626,99633,99634)
<pre>
439463330 reads; of these:
  439463330 (100.00%) were unpaired; of these:
    5393540 (1.23%) aligned 0 times
    328159166 (74.67%) aligned exactly 1 time
    105910624 (24.10%) aligned >1 times
98.77% overall alignment rate
[bam_sort_core] merging from 139 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 3:33:05. Running peak memory: 4.058GB.  
  PID: 99626;	Command: bowtie2;	Return code: 0;	Memory used: 4.058GB  
  PID: 99633;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 99634;	Command: samtools;	Return code: 0;	Memory used: 0.962GB


> `samtools view -q 10 -b -@ 24 -U /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam` (127236)
<pre>
</pre>
Command completed. Elapsed time: 0:09:45. Running peak memory: 4.058GB.  
  PID: 127236;	Command: samtools;	Return code: 0;	Memory used: 0.029GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	434069790	PEPPRO	_RES_

> `QC_filtered_reads`	46741044	PEPPRO	_RES_

> `Aligned_reads`	387328746	PEPPRO	_RES_

> `Alignment_rate`	80.0	PEPPRO	_RES_

> `Total_efficiency`	78.0	PEPPRO	_RES_

> `Read_depth`	29.61	PEPPRO	_RES_

### Compress all unmapped read files (05-27 15:30:17) elapsed: 16838.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/prealignments/K562_PRO-seq_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 24 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/prealignments/K562_PRO-seq_human_rDNA_unmap.fq` (134857)
<pre>
</pre>
Command completed. Elapsed time: 0:07:51. Running peak memory: 4.058GB.  
  PID: 134857;	Command: pigz;	Return code: 0;	Memory used: 0.019GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam` (135609)
<pre>
</pre>
Command completed. Elapsed time: 0:06:36. Running peak memory: 4.058GB.  
  PID: 135609;	Command: samtools;	Return code: 0;	Memory used: 0.027GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	9149071	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam` (136323)
<pre>
</pre>
Command completed. Elapsed time: 0:05:44. Running peak memory: 4.058GB.  
  PID: 136323;	Command: samtools;	Return code: 0;	Memory used: 0.024GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/chr_sizes.bed` (137408,137409,137410,137411)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.058GB.  
  PID: 137410;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 137408;	Command: samtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 137411;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 137409;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `samtools view -L /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/chr_sizes.bed -b -@ 24 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_noMT.bam` (137413)
<pre>
</pre>
Command completed. Elapsed time: 0:05:36. Running peak memory: 4.058GB.  
  PID: 137413;	Command: samtools;	Return code: 0;	Memory used: 0.03GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam` (137984)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 4.058GB.  
  PID: 137984;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam` (137993)
<pre>
</pre>
Command completed. Elapsed time: 0:05:37. Running peak memory: 4.058GB.  
  PID: 137993;	Command: samtools;	Return code: 0;	Memory used: 0.024GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (05-27 16:18:15) elapsed: 2878.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -c 24 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv` (141187)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/tmp_K562_PRO-seq_sort_my3ppij2'
Processing with 24 cores...
Discarding 81 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1']
Keeping 114 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270411v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:06:22. Running peak memory: 12.165GB.  
  PID: 141187;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 12.165GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_bamQC.tsv`

> `NRF`	0.28	PEPPRO	_RES_

> `PBC1`	0.54	PEPPRO	_RES_

> `PBC2`	3.89	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_unmap.bam`  

> `samtools view -b -@ 24 -f 4  /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_unmap.bam` (142117)
<pre>
</pre>
Command completed. Elapsed time: 0:01:19. Running peak memory: 12.165GB.  
  PID: 142117;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools view -c -f 4 -@ 24 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_temp.bam`

> `Unmapped_reads`	5393540	PEPPRO	_RES_

### Split BAM by strand (05-27 16:27:04) elapsed: 529.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam` (142526)
<pre>
</pre>
Command completed. Elapsed time: 0:28:09. Running peak memory: 12.165GB.  
  PID: 142526;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam` (146543)
<pre>
</pre>
Command completed. Elapsed time: 0:26:28. Running peak memory: 12.165GB.  
  PID: 146543;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (05-27 17:21:42) elapsed: 3278.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (150242)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.165GB.  
  PID: 150242;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/plus_TSS.tsv -p ends -c 24 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_plus_TssEnrichment.txt` (150249)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 12.165GB.  
  PID: 150249;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.25GB


> `TSS_coding_score`	13.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/minus_TSS.tsv -p ends -c 24 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_minus_TssEnrichment.txt` (150345)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 12.165GB.  
  PID: 150345;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.338GB


> `TSS_non-coding_score`	5.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_minus_TssEnrichment.txt` (150424)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 12.165GB.  
  PID: 150424;	Command: Rscript;	Return code: 0;	Memory used: 0.203GB

> `TSS enrichment`	QC_hg38/K562_PRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt` (150455,150456,150457,150458)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.165GB.  
  PID: 150455;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 150457;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 150456;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 150458;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt` (150460)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.165GB.  
  PID: 150460;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (05-27 17:23:02) elapsed: 81.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_tss.bed` (150462,150463)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 12.165GB.  
  PID: 150462;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 150463;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed` (150469,150470)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.165GB.  
  PID: 150469;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 150470;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_TSS_density.bed` (150472,150473,150474,150475)
<pre>
</pre>
Command completed. Elapsed time: 0:08:25. Running peak memory: 12.165GB.  
  PID: 150473;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 150475;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 150472;	Command: bedtools;	Return code: 0;	Memory used: 0.103GB  
  PID: 150474;	Command: sort;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_gene_body_density.bed` (151533,151538,151539)
<pre>
</pre>
Command completed. Elapsed time: 0:14:05. Running peak memory: 12.165GB.  
  PID: 151538;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 151533;	Command: bedtools;	Return code: 0;	Memory used: 1.018GB  
  PID: 151539;	Command: sort;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed` (153123,153129,153130)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.165GB.  
  PID: 153123;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 153130;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 153129;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed` (153135)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 12.165GB.  
  PID: 153135;	Command: Rscript;	Return code: 0;	Memory used: 0.32GB

> `Pause index`	QC_hg38/K562_PRO-seq_pause_index.pdf	Pause index	QC_hg38/K562_PRO-seq_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed.gz`  

> `pigz -f -p 24 -f /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_pause_index.bed` (153162)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.165GB.  
  PID: 153162;	Command: pigz;	Return code: 0;	Memory used: 0.004GB


### Calculate Fraction of Reads in pre-mature mRNA (05-27 17:45:39) elapsed: 1357.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam`
387328746 138157932

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam`
387328746 131399421

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_gene_sort.bed` (154423,154428)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.165GB.  
  PID: 154423;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 154428;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_gene_coverage.bed` (154430)
<pre>
</pre>
Command completed. Elapsed time: 0:14:24. Running peak memory: 12.165GB.  
  PID: 154430;	Command: bedtools;	Return code: 0;	Memory used: 1.018GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed.gz` (156126)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.165GB.  
  PID: 156126;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 24 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed` (156132)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.165GB.  
  PID: 156132;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (05-27 18:11:15) elapsed: 1536.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/raw/hg38_annotations.bed` (156145)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.165GB.  
  PID: 156145;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer_sort.bed` (156148,156149,156150,156151)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.165GB.  
  PID: 156148;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 156149;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 156151;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 156150;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_plus_coverage.bed` (156153)
<pre>
</pre>
Command completed. Elapsed time: 0:04:24. Running peak memory: 12.165GB.  
  PID: 156153;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_minus_coverage.bed` (156618)
<pre>
</pre>
Command completed. Elapsed time: 0:04:18. Running peak memory: 12.165GB.  
  PID: 156618;	Command: bedtools;	Return code: 0;	Memory used: 0.068GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_sort.bed` (156881,156882,156883,156884)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.165GB.  
  PID: 156881;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 156882;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 156884;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 156883;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_plus_coverage.bed` (156886)
<pre>
</pre>
Command completed. Elapsed time: 0:04:46. Running peak memory: 12.165GB.  
  PID: 156886;	Command: bedtools;	Return code: 0;	Memory used: 0.59GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_minus_coverage.bed` (157445)
<pre>
</pre>
Command completed. Elapsed time: 0:04:38. Running peak memory: 12.165GB.  
  PID: 157445;	Command: bedtools;	Return code: 0;	Memory used: 0.471GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region"` (157933)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.165GB.  
  PID: 157933;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed` (157934,157935,157936,157937)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.165GB.  
  PID: 157934;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 157936;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 157935;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 157937;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (157940)
<pre>
</pre>
Command completed. Elapsed time: 0:04:10. Running peak memory: 12.165GB.  
  PID: 157940;	Command: bedtools;	Return code: 0;	Memory used: 0.098GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_minus_coverage.bed` (158514)
<pre>
</pre>
Command completed. Elapsed time: 0:04:04. Running peak memory: 12.165GB.  
  PID: 158514;	Command: bedtools;	Return code: 0;	Memory used: 0.165GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR"` (158962)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.165GB.  
  PID: 158962;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR_sort.bed` (158963,158964,158965,158966)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.165GB.  
  PID: 158963;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 158964;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 158966;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 158965;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_plus_coverage.bed` (158969)
<pre>
</pre>
Command completed. Elapsed time: 0:04:03. Running peak memory: 12.165GB.  
  PID: 158969;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_minus_coverage.bed` (159485)
<pre>
</pre>
Command completed. Elapsed time: 0:03:56. Running peak memory: 12.165GB.  
  PID: 159485;	Command: bedtools;	Return code: 0;	Memory used: 0.064GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR"` (159929)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.165GB.  
  PID: 159929;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR_sort.bed` (159930,159931,159932,159933)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.165GB.  
  PID: 159930;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 159931;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 159933;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 159932;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_plus_coverage.bed` (159936)
<pre>
</pre>
Command completed. Elapsed time: 0:04:09. Running peak memory: 12.165GB.  
  PID: 159936;	Command: bedtools;	Return code: 0;	Memory used: 0.104GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_minus_coverage.bed` (160185)
<pre>
</pre>
Command completed. Elapsed time: 0:04:03. Running peak memory: 12.165GB.  
  PID: 160185;	Command: bedtools;	Return code: 0;	Memory used: 0.162GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Exon_sort.bed` (160697,160698,160699,160700)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 12.165GB.  
  PID: 160697;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 160698;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 160700;	Command: bedtools;	Return code: 0;	Memory used: 0.179GB  
  PID: 160699;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_plus_coverage.bed` (160707)
<pre>
</pre>
Command completed. Elapsed time: 0:04:30. Running peak memory: 12.165GB.  
  PID: 160707;	Command: bedtools;	Return code: 0;	Memory used: 0.465GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_minus_coverage.bed` (161187)
<pre>
</pre>
Command completed. Elapsed time: 0:04:25. Running peak memory: 12.165GB.  
  PID: 161187;	Command: bedtools;	Return code: 0;	Memory used: 0.185GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Intron_sort.bed` (161780,161781,161782,161783)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.165GB.  
  PID: 161780;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 161782;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 161781;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 161783;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_plus_coverage.bed` (161787)
<pre>
</pre>
Command completed. Elapsed time: 0:05:31. Running peak memory: 12.165GB.  
  PID: 161787;	Command: bedtools;	Return code: 0;	Memory used: 0.451GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_minus_coverage.bed` (162332)
<pre>
</pre>
Command completed. Elapsed time: 0:05:36. Running peak memory: 12.165GB.  
  PID: 162332;	Command: bedtools;	Return code: 0;	Memory used: 0.431GB


### Plot cFRiF/FRiF (05-27 19:13:58) elapsed: 3763.0 _TIME_


> `samtools view -@ 24 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq -z 3099922541 -n 192750289 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_plus_coverage.bed` (163007)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 12.165GB.  
  PID: 163007;	Command: Rscript;	Return code: 0;	Memory used: 0.436GB

> `cFRiF`	QC_hg38/K562_PRO-seq_cFRiF.pdf	cFRiF	QC_hg38/K562_PRO-seq_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq -z 3099922541 -n 192750289 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_Intron_plus_coverage.bed` (163254)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 12.165GB.  
  PID: 163254;	Command: Rscript;	Return code: 0;	Memory used: 0.492GB

> `FRiF`	QC_hg38/K562_PRO-seq_FRiF.pdf	FRiF	QC_hg38/K562_PRO-seq_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (05-27 19:15:33) elapsed: 95.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_exons_sort.bed` (163299,163301)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 12.165GB.  
  PID: 163301;	Command: bedtools;	Return code: 0;	Memory used: 0.097GB  
  PID: 163299;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_introns_sort.bed` (163306,163308,163309)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 12.165GB.  
  PID: 163306;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 163309;	Command: bedtools;	Return code: 0;	Memory used: 0.003GB  
  PID: 163308;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_coverage.bed` (163315)
<pre>
</pre>
Command completed. Elapsed time: 0:09:24. Running peak memory: 12.165GB.  
  PID: 163315;	Command: bedtools;	Return code: 0;	Memory used: 0.291GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_coverage.bed` (164354)
<pre>
</pre>
Command completed. Elapsed time: 0:12:04. Running peak memory: 12.165GB.  
  PID: 164354;	Command: bedtools;	Return code: 0;	Memory used: 0.485GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/387.328746)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_rpkm.bed` (165623,165630,165631)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.165GB.  
  PID: 165623;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 165631;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 165630;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/387.328746)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_rpkm.bed` (165634,165635,165636)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 12.165GB.  
  PID: 165634;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 165636;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 165635;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed` (165639,165640,165641)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.165GB.  
  PID: 165639;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 165641;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 165640;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed --annotate` (165647)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 12.165GB.  
  PID: 165647;	Command: Rscript;	Return code: 0;	Memory used: 0.32GB

> `mRNA contamination`	QC_hg38/K562_PRO-seq_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_PRO-seq_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed.gz`  

> `pigz -f -p 24 -f /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/QC_hg38/K562_PRO-seq_exon_intron_ratios.bed` (165674)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 12.165GB.  
  PID: 165674;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (05-27 19:37:17) elapsed: 1305.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam` (165684)
<pre>
</pre>
Command completed. Elapsed time: 0:02:53. Running peak memory: 12.165GB.  
  PID: 165684;	Command: samtools;	Return code: 0;	Memory used: 0.02GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_smooth_body_0-mer.bw -p 16 --variable-step --tail-edge` (166147)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_plus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_plus_cuttrace_ccjcwirk'
Processing with 8 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 90 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1']
Keeping 105 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 105 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_exact_body_0-mer.bw'
Merging 105 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:10:13. Running peak memory: 12.165GB.  
  PID: 166147;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 4.32GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam` (168231)
<pre>
</pre>
Command completed. Elapsed time: 0:02:47. Running peak memory: 12.165GB.  
  PID: 168231;	Command: samtools;	Return code: 0;	Memory used: 0.021GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_smooth_body_0-mer.bw -p 16 --variable-step --tail-edge` (168401)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/aligned_hg38/K562_PRO-seq_minus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_minus_cuttrace_cko3imu0'
Processing with 8 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
stdin is empty of data
psutil.NoSuchProcess process no longer exists (pid=169063)
Warning: couldn't add memory use for process: 168401
stdin is empty of data
Discarding 93 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 102 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270411v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270336v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 102 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_exact_body_0-mer.bw'
Merging 102 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/K562_PRO-seq/signal_hg38/K562_PRO-seq_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:52. Running peak memory: 12.165GB.  
  PID: 168401;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.311GB

Starting cleanup: 53 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  11:48:14
*  Total elapsed time (all runs):  15:14:24
*         Peak memory (this run):  12.1652 GB
*        Pipeline completed time: 2020-05-27 20:03:12
