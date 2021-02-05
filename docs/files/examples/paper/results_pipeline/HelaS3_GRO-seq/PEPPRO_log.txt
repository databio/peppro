### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name HelaS3_GRO-seq --genome hg38 --input /project/shefflab/data//sra_fastq/SRR1693611.fastq.gz /project/shefflab/data//sra_fastq/SRR1693612.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol GRO --umi-len 0 --scale --prealignments human_rDNA --dirty`
*         Compute host:  udc-aj37-18c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/
*  Pipeline started at:   (06-11 17:32:09) elapsed: 11.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//sra_fastq/SRR1693611.fastq.gz', '/project/shefflab/data//sra_fastq/SRR1693612.fastq.gz']`
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
*        `sample_name`:  `HelaS3_GRO-seq`
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

Local input file: /project/shefflab/data//sra_fastq/SRR1693611.fastq.gz

> `File_mb`	3550.33	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected GRO input

### Merge/link and fastq conversion:  (06-11 17:32:11) elapsed: 2.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz`  

> `cat /project/shefflab/data//sra_fastq/SRR1693611.fastq.gz /project/shefflab/data//sra_fastq/SRR1693612.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz` (51527)
<pre>
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 0.003GB.  
  PID: 51527;	Command: cat;	Return code: 0;	Memory used: 0.003GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1.fastq` (51577)
<pre>
</pre>
Command completed. Elapsed time: 0:02:15. Running peak memory: 0.003GB.  
  PID: 51577;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	72718941	PEPPRO	_RES_

> `Fastq_reads`	72718941	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz']

### FASTQ processing:  (06-11 17:38:35) elapsed: 384.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_processed.fastq`  

> `(cutadapt -j 12 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt` (52260)
<pre>
</pre>
Command completed. Elapsed time: 0:03:25. Running peak memory: 4.586GB.  
  PID: 52260;	Command: cutadapt;	Return code: 0;	Memory used: 4.586GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_noadap.fastq | seqtk seq -L 2 - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_processed.fastq` (52745,52746)
<pre>
</pre>
Command completed. Elapsed time: 0:02:38. Running peak memory: 4.586GB.  
  PID: 52745;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 52746;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	22179715.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	0.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	0.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastqc/HelaS3_GRO-seq_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (54520)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.586GB.  
  PID: 54520;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	72718941	PEPPRO	_RES_

> `Trim_loss_rate`	0.0	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_processed.fastq` (55987)
<pre>
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
Command completed. Elapsed time: 0:03:57. Running peak memory: 4.586GB.  
  PID: 55987;	Command: fastqc;	Return code: 0;	Memory used: 0.192GB

> `FastQC report r1`	fastqc/HelaS3_GRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/processed_R1.flag` (57324)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.586GB.  
  PID: 57324;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Plot adapter insertion distribution (06-11 17:54:51) elapsed: 976.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt` (57326)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 4.586GB.  
  PID: 57326;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Adapter insertion distribution`	cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	0	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 17:54:58) elapsed: 7.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk 'END {print $1}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 10 && $1 >= 14){degradedSum += $2}; ($1 >= 4 && $1 <= 0){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.0	PEPPRO	_RES_

### Prealignments (06-11 17:54:58) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 17:54:58) elapsed: 0.0 _TIME_


> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id HelaS3_GRO-seq -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
72718941 reads; of these:
  72718941 (100.00%) were unpaired; of these:
    31839185 (43.78%) aligned 0 times
    40879756 (56.22%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
56.22% overall alignment rate

> `Aligned_reads_human_rDNA`	40879756.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	56.22	PEPPRO	_RES_

### Map to genome (06-11 18:03:13) elapsed: 495.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam`  

> `bowtie2 -p 12 --very-sensitive --rg-id HelaS3_GRO-seq -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/tmplwu3sau6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam` (59918,59923,59924)
<pre>
31839185 reads; of these:
  31839185 (100.00%) were unpaired; of these:
    10809494 (33.95%) aligned 0 times
    13023651 (40.90%) aligned exactly 1 time
    8006040 (25.15%) aligned >1 times
66.05% overall alignment rate
[bam_sort_core] merging from 8 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:12:23. Running peak memory: 4.586GB.  
  PID: 59918;	Command: bowtie2;	Return code: 0;	Memory used: 3.666GB  
  PID: 59923;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 59924;	Command: samtools;	Return code: 0;	Memory used: 0.899GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam` (63877)
<pre>
</pre>
Command completed. Elapsed time: 0:01:41. Running peak memory: 4.586GB.  
  PID: 63877;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	21029691	PEPPRO	_RES_

> `QC_filtered_reads`	6163527	PEPPRO	_RES_

> `Aligned_reads`	14866164	PEPPRO	_RES_

> `Alignment_rate`	20.44	PEPPRO	_RES_

> `Total_efficiency`	20.44	PEPPRO	_RES_

> `Read_depth`	2.34	PEPPRO	_RES_

### Compress all unmapped read files (06-11 18:26:53) elapsed: 1420.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq` (67152)
<pre>
</pre>
Command completed. Elapsed time: 0:01:03. Running peak memory: 4.586GB.  
  PID: 67152;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam` (67418)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 4.586GB.  
  PID: 67418;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	1906316	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam` (67506)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 4.586GB.  
  PID: 67506;	Command: samtools;	Return code: 0;	Memory used: 0.007GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/chr_sizes.bed` (67553,67554,67555,67556)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.586GB.  
  PID: 67553;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 67555;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 67554;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 67556;	Command: grep;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_noMT.bam` (67558)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 4.586GB.  
  PID: 67558;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam` (67633)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.586GB.  
  PID: 67633;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam` (67637)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 4.586GB.  
  PID: 67637;	Command: samtools;	Return code: 0;	Memory used: 0.007GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	51	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-11 18:29:41) elapsed: 168.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv` (67810)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/tmp_HelaS3_GRO-seq_sort_whhdf766'
Processing with 12 cores...
Discarding 117 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr22_KI270735v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 78 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270538v1', 'chrUn_KI270582v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000224v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 4.586GB.  
  PID: 67810;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.267GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv`

> `NRF`	0.72	PEPPRO	_RES_

> `PBC1`	0.91	PEPPRO	_RES_

> `PBC2`	13.55	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_unmap.bam`  

> `samtools view -b -@ 12 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_unmap.bam` (68140)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 4.586GB.  
  PID: 68140;	Command: samtools;	Return code: 0;	Memory used: 0.017GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_temp.bam`

> `Unmapped_reads`	10809494	PEPPRO	_RES_

### Split BAM by strand (06-11 18:30:19) elapsed: 38.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam` (68322)
<pre>
</pre>
Command completed. Elapsed time: 0:00:57. Running peak memory: 4.586GB.  
  PID: 68322;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam` (68549)
<pre>
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 4.586GB.  
  PID: 68549;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-11 18:32:10) elapsed: 111.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (68757)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.586GB.  
  PID: 68757;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_plus_TssEnrichment.txt` (68761)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 4.586GB.  
  PID: 68761;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.509GB


> `TSS_coding_score`	7.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_minus_TssEnrichment.txt` (68819)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.586GB.  
  PID: 68819;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.863GB


> `TSS_non-coding_score`	4.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_minus_TssEnrichment.txt` (68863)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 4.586GB.  
  PID: 68863;	Command: Rscript;	Return code: 0;	Memory used: 0.084GB

> `TSS enrichment`	QC_hg38/HelaS3_GRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/HelaS3_GRO-seq_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt` (68914,68915,68916,68917)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.586GB.  
  PID: 68914;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 68916;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 68915;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 68917;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt` (68922)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.586GB.  
  PID: 68922;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-11 18:32:36) elapsed: 26.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_tss.bed` (68924,68925)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 4.586GB.  
  PID: 68924;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 68925;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed` (68932,68933)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.586GB.  
  PID: 68932;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 68933;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSS_density.bed` (68936,68937,68938,68939)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 4.586GB.  
  PID: 68936;	Command: bedtools;	Return code: 0;	Memory used: 0.056GB  
  PID: 68938;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 68937;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 68939;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_gene_body_density.bed` (68995,68996,68997)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 4.586GB.  
  PID: 68996;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 68995;	Command: bedtools;	Return code: 0;	Memory used: 0.137GB  
  PID: 68997;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed`  
awk: cmd. line:1: { if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}
awk: cmd. line:1:             ^ syntax error
awk: cmd. line:1: { if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}
awk: cmd. line:1:                                                                                                                      ^ syntax error

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/tmp0l2qc6_9` (69034,69035,69036)
<pre>
</pre>
Got SIGTERM. Failing gracefully... (06-11 18:55:15) elapsed: 1359.0 _TIME_
Child process 69034 (join) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/recover.lock.QC_hg38__HelaS3_GRO-seq_pause_index.bed

### Pipeline failed at:  (06-11 18:55:15) elapsed: 0.0 _TIME_

Total time: 1:23:16
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name HelaS3_GRO-seq --genome hg38 --input /project/shefflab/data//sra_fastq/SRR1693611.fastq.gz /project/shefflab/data//sra_fastq/SRR1693612.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol GRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba27-12
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/
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
*              `input`:  `['/project/shefflab/data//sra_fastq/SRR1693611.fastq.gz', '/project/shefflab/data//sra_fastq/SRR1693612.fastq.gz']`
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
*        `sample_name`:  `HelaS3_GRO-seq`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//sra_fastq/SRR1693611.fastq.gz

> `File_mb`	3550.33	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected GRO input

### Merge/link and fastq conversion:  (06-11 19:10:15) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz']

### FASTQ processing:  (06-11 19:10:15) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-11 19:10:15) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-11 19:10:15) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 19:10:15) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq

### Map to genome (06-11 19:10:15) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam`  

### Compress all unmapped read files (06-11 19:10:15) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-11 19:10:15) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_unmap.bam`  

### Split BAM by strand (06-11 19:10:15) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam`  

### Calculate TSS enrichment (06-11 19:10:15) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/HelaS3_GRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/HelaS3_GRO-seq_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-11 19:10:15) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_gene_body_density.bed`  
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/lock.QC_hg38__HelaS3_GRO-seq_pause_index.bed
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/tmpv8eq039_` (115692,115693,115694)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.006GB.  
  PID: 115692;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 115694;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 115693;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/tmpv8eq039_ | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-11 19:10:16) elapsed: 0.0 _TIME_

Total time: 0:00:01
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/tmpv8eq039_ | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name HelaS3_GRO-seq --genome hg38 --input /project/shefflab/data//sra_fastq/SRR1693611.fastq.gz /project/shefflab/data//sra_fastq/SRR1693612.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol GRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aj37-18c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/
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
*              `input`:  `['/project/shefflab/data//sra_fastq/SRR1693611.fastq.gz', '/project/shefflab/data//sra_fastq/SRR1693612.fastq.gz']`
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
*        `sample_name`:  `HelaS3_GRO-seq`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//sra_fastq/SRR1693611.fastq.gz

> `File_mb`	3550.33	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected GRO input

### Merge/link and fastq conversion:  (06-14 21:14:11) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/HelaS3_GRO-seq_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/HelaS3_GRO-seq.merged.fastq.gz']

### FASTQ processing:  (06-14 21:14:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-14 21:14:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/HelaS3_GRO-seq_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-14 21:14:11) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:14:11) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq

### Map to genome (06-14 21:14:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam`  

### Compress all unmapped read files (06-14 21:14:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/prealignments/HelaS3_GRO-seq_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:14:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_unmap.bam`  

### Split BAM by strand (06-14 21:14:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam`  

### Calculate TSS enrichment (06-14 21:14:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/HelaS3_GRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/HelaS3_GRO-seq_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:14:11) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/tmpdxowb6b2` (428182,428190,428202)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.005GB.  
  PID: 428182;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 428202;	Command: env;	Return code: 0;	Memory used: 0.005GB  
  PID: 428190;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/tmpdxowb6b2 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0143043) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/tmpdxowb6b2 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed` (428270)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.005GB.  
  PID: 428270;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	4.92	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed` (428309)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.279GB.  
  PID: 428309;	Command: Rscript;	Return code: 0;	Memory used: 0.279GB

> `Pause index`	QC_hg38/HelaS3_GRO-seq_pause_index.pdf	Pause index	QC_hg38/HelaS3_GRO-seq_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_pause_index.bed` (429785)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.279GB.  
  PID: 429785;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:14:18) elapsed: 6.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam`
14866164 5457524

> `Plus_FRiP`	0.37	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam`
14866164 4922927

> `Minus_FRiP`	0.33	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_gene_sort.bed` (430831,430832)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.279GB.  
  PID: 430831;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 430832;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_gene_coverage.bed` (430834)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 0.279GB.  
  PID: 430834;	Command: bedtools;	Return code: 0;	Memory used: 0.104GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed.gz` (431057)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.279GB.  
  PID: 431057;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed` (431058)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.279GB.  
  PID: 431058;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:15:07) elapsed: 49.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/raw/hg38_annotations.bed` (431067)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.279GB.  
  PID: 431067;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer_sort.bed` (431070,431071,431072,431073)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.279GB.  
  PID: 431070;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 431071;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 431073;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB  
  PID: 431072;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_plus_coverage.bed` (431075)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 0.279GB.  
  PID: 431075;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_minus_coverage.bed` (431110)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 0.279GB.  
  PID: 431110;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_sort.bed` (431123,431124,431125,431126)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.279GB.  
  PID: 431123;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 431124;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 431126;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 431125;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_plus_coverage.bed` (431128)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 0.279GB.  
  PID: 431128;	Command: bedtools;	Return code: 0;	Memory used: 0.061GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_minus_coverage.bed` (431141)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 0.279GB.  
  PID: 431141;	Command: bedtools;	Return code: 0;	Memory used: 0.067GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region"` (431156)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.279GB.  
  PID: 431156;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed` (431157,431158,431159,431160)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.279GB.  
  PID: 431157;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 431159;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 431158;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 431160;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (431163)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 0.279GB.  
  PID: 431163;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_minus_coverage.bed` (431180)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 0.279GB.  
  PID: 431180;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR"` (431212)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.279GB.  
  PID: 431212;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR_sort.bed` (431213,431214,431215,431216)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.279GB.  
  PID: 431213;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 431214;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 431216;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 431215;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_plus_coverage.bed` (431219)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 0.279GB.  
  PID: 431219;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_minus_coverage.bed` (431231)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 0.279GB.  
  PID: 431231;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR"` (431240)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.279GB.  
  PID: 431240;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR_sort.bed` (431241,431242,431243,431244)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.279GB.  
  PID: 431241;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 431242;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 431244;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 431243;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_plus_coverage.bed` (431246)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 0.279GB.  
  PID: 431246;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_minus_coverage.bed` (431264)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 0.279GB.  
  PID: 431264;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon_sort.bed` (431275,431276,431277,431278)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.279GB.  
  PID: 431275;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 431276;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 431278;	Command: bedtools;	Return code: 0;	Memory used: 0.159GB  
  PID: 431277;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_plus_coverage.bed` (431282)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 0.279GB.  
  PID: 431282;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_minus_coverage.bed` (431333)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 0.279GB.  
  PID: 431333;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron_sort.bed` (436224,436225,436226,436227)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.279GB.  
  PID: 436224;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 436226;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 436225;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 436227;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_plus_coverage.bed` (437383)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 0.279GB.  
  PID: 437383;	Command: bedtools;	Return code: 0;	Memory used: 0.04GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_minus_coverage.bed` (453287)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 0.279GB.  
  PID: 453287;	Command: bedtools;	Return code: 0;	Memory used: 0.05GB


### Plot cFRiF/FRiF (06-14 21:17:29) elapsed: 142.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s HelaS3_GRO-seq -z 3099922541 -n 7116450 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_plus_coverage.bed` (11157)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 0.475GB.  
  PID: 11157;	Command: Rscript;	Return code: 0;	Memory used: 0.475GB

> `cFRiF`	QC_hg38/HelaS3_GRO-seq_cFRiF.pdf	cFRiF	QC_hg38/HelaS3_GRO-seq_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s HelaS3_GRO-seq -z 3099922541 -n 7116450 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_Intron_plus_coverage.bed` (22651)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 0.513GB.  
  PID: 22651;	Command: Rscript;	Return code: 0;	Memory used: 0.513GB

> `FRiF`	QC_hg38/HelaS3_GRO-seq_FRiF.pdf	FRiF	QC_hg38/HelaS3_GRO-seq_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 21:18:27) elapsed: 58.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_exons_sort.bed` (44531,44533)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.513GB.  
  PID: 44533;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB  
  PID: 44531;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_introns_sort.bed` (50607,50620,50621)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.513GB.  
  PID: 50607;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 50621;	Command: bedtools;	Return code: 0;	Memory used: 0.003GB  
  PID: 50620;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_coverage.bed` (55318)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 0.513GB.  
  PID: 55318;	Command: bedtools;	Return code: 0;	Memory used: 0.397GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_coverage.bed` (61248)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 0.513GB.  
  PID: 61248;	Command: bedtools;	Return code: 0;	Memory used: 0.091GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/14.866164)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_rpkm.bed` (61372,61373,61374)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.513GB.  
  PID: 61372;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 61374;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 61373;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/14.866164)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_rpkm.bed` (61376,61377,61378)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.513GB.  
  PID: 61376;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 61378;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 61377;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed` (61381,61382,61383)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.513GB.  
  PID: 61381;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 61383;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 61382;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	3.22	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed --annotate` (61389)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.513GB.  
  PID: 61389;	Command: Rscript;	Return code: 0;	Memory used: 0.317GB

> `mRNA contamination`	QC_hg38/HelaS3_GRO-seq_mRNA_contamination.pdf	mRNA contamination	QC_hg38/HelaS3_GRO-seq_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/QC_hg38/HelaS3_GRO-seq_exon_intron_ratios.bed` (61409)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.513GB.  
  PID: 61409;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (06-14 21:19:21) elapsed: 54.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam` (61417)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.513GB.  
  PID: 61417;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_smooth_body_0-mer.bw -p 8 --variable-step --scale 14866164.0` (61423)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_plus.bam'
Temporary files will be stored in: 'tmp_HelaS3_GRO-seq_plus_cuttrace_dqxcnk7g'
Processing with 4 cores...
Discarding 134 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr14_GL000225v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 61 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 61 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_exact_body_0-mer.bw'
Merging 61 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:05:23. Running peak memory: 3.503GB.  
  PID: 61423;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.503GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam` (151408)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.503GB.  
  PID: 151408;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_smooth_body_0-mer.bw -p 8 --variable-step --scale 14866164.0` (151418)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/aligned_hg38/HelaS3_GRO-seq_minus.bam'
Temporary files will be stored in: 'tmp_HelaS3_GRO-seq_minus_cuttrace_a8m5gu0r'
Processing with 4 cores...
stdin is empty of data
Discarding 135 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_GL000216v2', 'chrEBV']
Keeping 60 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270538v1', 'chrUn_KI270582v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270742v1', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 60 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_exact_body_0-mer.bw'
Merging 60 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HelaS3_GRO-seq/signal_hg38/HelaS3_GRO-seq_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:05:18. Running peak memory: 3.524GB.  
  PID: 151418;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.524GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:16:11
*  Total elapsed time (all runs):  1:08:54
*         Peak memory (this run):  3.524 GB
*        Pipeline completed time: 2020-06-14 21:30:12
