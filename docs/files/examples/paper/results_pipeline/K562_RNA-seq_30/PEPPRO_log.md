### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name K562_RNA-seq_30 --genome hg38 --input /project/shefflab/data/guertin/fastq/K562_30pct_RNArc_r2.fastq.gz --single-or-paired single --protocol PRO --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/dev4/results_pipeline -P 4 -M 16000`
*         Compute host:  udc-ba27-32c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/
*  Pipeline started at:   (02-27 09:23:37) elapsed: 1.0 _TIME_

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
*              `cores`:  `4`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `False`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data/guertin/fastq/K562_30pct_RNArc_r2.fastq.gz']`
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
*           `protocol`:  `PRO`
*            `recover`:  `False`
*        `sample_name`:  `K562_RNA-seq_30`
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

Local input file: /project/shefflab/data/guertin/fastq/K562_30pct_RNArc_r2.fastq.gz

> `File_mb`	803.62	PEPPRO	_RES_

> `Read_type`	single	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-27 09:23:37) elapsed: 0.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/K562_30pct_RNArc_r2.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz` (86144)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.002GB.  
  PID: 86144;	Command: ln;	Return code: 0;	Memory used: 0.002GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1.fastq`  

> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1.fastq` (86145)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 0.003GB.  
  PID: 86145;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	10000000	PEPPRO	_RES_

> `Fastq_reads`	10000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/raw/K562_RNA-seq_30.fastq.gz']

### FASTQ processing:  (02-27 09:24:08) elapsed: 31.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_processed.fastq`  

> `(cutadapt -j 4 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt` (86204)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 0.172GB.  
  PID: 86204;	Command: cutadapt;	Return code: 0;	Memory used: 0.172GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_processed.fastq` (86238,86239)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 0.172GB.  
  PID: 86238;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 86239;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	7013383.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Uninformative_adapter_reads`	174726.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	1.7473	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/fastqc/K562_RNA-seq_30_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (86273)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.172GB.  
  PID: 86273;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	9825274	PEPPRO	_RES_

> `Trim_loss_rate`	1.75	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_processed.fastq` (86279)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
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
Command completed. Elapsed time: 0:00:36. Running peak memory: 0.192GB.  
  PID: 86279;	Command: fastqc;	Return code: 0;	Memory used: 0.192GB

> `FastQC report r1`	fastqc/K562_RNA-seq_30_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/fastq/processed_R1.flag` (86521)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.192GB.  
  PID: 86521;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (02-27 09:25:27) elapsed: 79.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/cutadapt` (86522)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.201GB.  
  PID: 86522;	Command: Rscript;	Return code: 0;	Memory used: 0.201GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_30_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-27 09:25:33) elapsed: 6.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/cutadapt/K562_RNA-seq_30_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2303	PEPPRO	_RES_

### Prealignments (02-27 09:25:33) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-27 09:25:33) elapsed: 0.0 _TIME_


> `(bowtie2 -p 4 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_30 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/fastq/K562_RNA-seq_30_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
9825274 reads; of these:
  9825274 (100.00%) were unpaired; of these:
    9122083 (92.84%) aligned 0 times
    703191 (7.16%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
7.16% overall alignment rate

> `Aligned_reads_human_rDNA`	703191.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	7.16	PEPPRO	_RES_

### Map to genome (02-27 09:26:44) elapsed: 70.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam`  

> `bowtie2 -p 4 --very-sensitive -X 2000 --rg-id K562_RNA-seq_30 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/tmp1yhq07vf -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam` (86638,86639,86640)
<pre>
9122083 reads; of these:
  9122083 (100.00%) were unpaired; of these:
    524877 (5.75%) aligned 0 times
    6088436 (66.74%) aligned exactly 1 time
    2508770 (27.50%) aligned >1 times
94.25% overall alignment rate
[bam_sort_core] merging from 2 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:11:14. Running peak memory: 3.534GB.  
  PID: 86638;	Command: bowtie2;	Return code: 0;	Memory used: 3.534GB  
  PID: 86639;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 86640;	Command: samtools;	Return code: 0;	Memory used: 0.866GB


> `samtools view -q 10 -b -@ 4 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_fail_qc.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam` (87694)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 3.534GB.  
  PID: 87694;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	8597206	PEPPRO	_RES_

> `QC_filtered_reads`	1214446	PEPPRO	_RES_

> `Aligned_reads`	7382760	PEPPRO	_RES_

> `Alignment_rate`	75.14	PEPPRO	_RES_

> `Total_efficiency`	73.83	PEPPRO	_RES_

> `Read_depth`	2.15	PEPPRO	_RES_

### Compress all unmapped read files (02-27 09:45:24) elapsed: 1120.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/prealignments/K562_RNA-seq_30_human_rDNA_unmap.fq` (88511)
<pre>
</pre>
Command completed. Elapsed time: 0:00:47. Running peak memory: 3.534GB.  
  PID: 88511;	Command: pigz;	Return code: 0;	Memory used: 0.004GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam` (88566)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.534GB.  
  PID: 88566;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	309258	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam` (88578)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 88578;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_noMT.bam` (88586,88587,88588,88589)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.534GB.  
  PID: 88588;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 88586;	Command: samtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 88587;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 88589;	Command: xargs;	Return code: 0;	Memory used: 0.029GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_noMT.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam` (88617)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 88617;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam` (88618)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 88618;	Command: samtools;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-27 09:47:17) elapsed: 113.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -c 4 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv` (88653)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/tmp_K562_RNA-seq_30_sort_2xfglb59'
Processing with 4 cores...
Discarding 115 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270731v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrEBV']
Keeping 80 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.534GB.  
  PID: 88653;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.435GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_bamQC.tsv`

> `NRF`	0.93	PEPPRO	_RES_

> `PBC1`	0.97	PEPPRO	_RES_

> `PBC2`	46.16	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_unmap.bam`  

> `samtools view -b -@ 4 -f 4  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_unmap.bam` (88681)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.534GB.  
  PID: 88681;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools view -c -f 4 -@ 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_temp.bam`

> `Unmapped_reads`	524877	PEPPRO	_RES_

### Split BAM by strand (02-27 09:47:42) elapsed: 25.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam` (88705)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 3.534GB.  
  PID: 88705;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam` (88735)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.534GB.  
  PID: 88735;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-27 09:48:52) elapsed: 70.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (88768)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 88768;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/plus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_plus_TssEnrichment.txt` (88769)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 88769;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.196GB


> `TSS_coding_score`	16.2	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/minus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_minus_TssEnrichment.txt` (88787)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 88787;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.2GB


> `TSS_non-coding_score`	5.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_minus_TssEnrichment.txt` (88803)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 88803;	Command: Rscript;	Return code: 0;	Memory used: 0.22GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_30_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_30_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt` (88820,88821,88822,88823)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 88820;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 88822;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 88821;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 88823;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt` (88825)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 88825;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-27 09:49:11) elapsed: 19.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_tss.bed` (88827,88828)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.534GB.  
  PID: 88827;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 88828;	Command: bedtools;	Return code: 0;	Memory used: 0.098GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_gene_body.bed` (88832,88833)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 88832;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 88833;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSS_density.bed` (88835,88836,88837,88838)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.534GB.  
  PID: 88837;	Command: sort;	Return code: 0;	Memory used: 0.012GB  
  PID: 88835;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 88838;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 88836;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_gene_body_density.bed` (88850,88851,88852)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.534GB.  
  PID: 88850;	Command: bedtools;	Return code: 0;	Memory used: 0.04GB  
  PID: 88852;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 88851;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_TSS_density.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed` (88864,88865,88866)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 88864;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 88866;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 88865;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	13.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed` (88871)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.534GB.  
  PID: 88871;	Command: Rscript;	Return code: 0;	Memory used: 0.369GB

> `Pause index`	QC_hg38/K562_RNA-seq_30_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_30_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_pause_index.bed` (88890)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 88890;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Calculate Fraction of Reads in pre-mature mRNA (02-27 09:49:42) elapsed: 32.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam`
7382760 2719315

> `Plus_FRiP`	0.37	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam`
7382760 2636063

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_gene_sort.bed` (88918,88919)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 88918;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 88919;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_gene_coverage.bed` (88922)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.534GB.  
  PID: 88922;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed.gz` (89146)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 89146;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed` (89147)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.534GB.  
  PID: 89147;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (02-27 09:50:11) elapsed: 28.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/raw/hg38_annotations.bed` (89164)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.534GB.  
  PID: 89164;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR"` (89166)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 89166;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR_sort.bed` (89167,89168,89169,89170)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 89167;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 89168;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 89170;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 89169;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_plus_coverage.bed` (89172)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 89172;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_minus_coverage.bed` (89179)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 89179;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR"` (89185)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 89185;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR_sort.bed` (89186,89187,89188,89189)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 89186;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 89187;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 89189;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 89188;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_plus_coverage.bed` (89192)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 89192;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_minus_coverage.bed` (89198)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 89198;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer_sort.bed` (89205,89206,89207,89208)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 89205;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 89206;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 89208;	Command: bedtools;	Return code: 0;	Memory used: 0.047GB  
  PID: 89207;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_plus_coverage.bed` (89210)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 89210;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_minus_coverage.bed` (89220)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 89220;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon_sort.bed` (89226,89227,89228,89229)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.534GB.  
  PID: 89226;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 89227;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 89229;	Command: bedtools;	Return code: 0;	Memory used: 0.174GB  
  PID: 89228;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_plus_coverage.bed` (89233)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 89233;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_minus_coverage.bed` (89240)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 89240;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron_sort.bed` (89248,89249,89250,89251)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 89248;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 89250;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 89249;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 89251;	Command: bedtools;	Return code: 0;	Memory used: 0.083GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_plus_coverage.bed` (89254)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 89254;	Command: bedtools;	Return code: 0;	Memory used: 0.023GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_minus_coverage.bed` (89261)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 89261;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_sort.bed` (89269,89270,89271,89272)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 89269;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 89270;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 89272;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 89271;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_plus_coverage.bed` (89274)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 89274;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_minus_coverage.bed` (89282)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 89282;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region"` (89288)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 89288;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed` (89289,89290,89291,89292)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 89289;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 89291;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 89290;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 89292;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_plus_coverage.bed` (89295)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 89295;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_minus_coverage.bed` (89302)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 89302;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB


### Plot cFRiF/FRiF (02-27 09:51:45) elapsed: 94.0 _TIME_


> `samtools view -@ 4 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_30 -z 3099922541 -n 3603001 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_plus_coverage.bed` (89318)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 3.534GB.  
  PID: 89318;	Command: Rscript;	Return code: 0;	Memory used: 0.477GB

> `cFRiF`	QC_hg38/K562_RNA-seq_30_cFRiF.pdf	cFRiF	QC_hg38/K562_RNA-seq_30_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_30 -z 3099922541 -n 3603001 -y frif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_Promoter_Flanking_Region_plus_coverage.bed` (89353)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 3.534GB.  
  PID: 89353;	Command: Rscript;	Return code: 0;	Memory used: 0.556GB

> `FRiF`	QC_hg38/K562_RNA-seq_30_FRiF.pdf	FRiF	QC_hg38/K562_RNA-seq_30_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-27 09:52:35) elapsed: 51.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_exons_sort.bed` (89382,89383)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 89383;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB  
  PID: 89382;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_introns_sort.bed` (89393,89394,89395)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 89393;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 89395;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 89394;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_coverage.bed` (89402)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.534GB.  
  PID: 89402;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_coverage.bed` (89413)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.534GB.  
  PID: 89413;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.38276)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_rpkm.bed` (89424,89425,89426)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 89424;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 89426;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 89425;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.38276)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_rpkm.bed` (89428,89429,89430)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 89428;	Command: awk;	Return code: 0;	Memory used: 0.006GB  
  PID: 89430;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 89429;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_introns_rpkm.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed` (89433,89434,89435)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 89433;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 89435;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 89434;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	3.21	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed --annotate` (89441)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.534GB.  
  PID: 89441;	Command: Rscript;	Return code: 0;	Memory used: 0.366GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_30_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_30_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/QC_hg38/K562_RNA-seq_30_exon_intron_ratios.bed` (89458)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 89458;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Produce bigWig files (02-27 09:53:15) elapsed: 40.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam` (89465)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.534GB.  
  PID: 89465;	Command: samtools;	Return code: 0;	Memory used: 0.005GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge` (89469)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_30_plus_cuttrace_cguw_64v'
Processing with 1 cores...
Discarding 126 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_GL000218v1', 'chrEBV']
Keeping 69 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270582v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2']
Reduce step (merge files)...
Merging 69 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_exact_body_0-mer.bw'
Merging 69 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:14:09. Running peak memory: 3.534GB.  
  PID: 89469;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.449GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam` (91511)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.534GB.  
  PID: 91511;	Command: samtools;	Return code: 0;	Memory used: 0.005GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge` (91515)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/aligned_hg38/K562_RNA-seq_30_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_30_minus_cuttrace_64mddjhd'
Processing with 1 cores...
Discarding 125 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr16_KI270728v1_random', 'chr22_KI270731v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 70 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 70 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_exact_body_0-mer.bw'
Merging 70 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_30/signal_hg38/K562_RNA-seq_30_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:13:54. Running peak memory: 3.534GB.  
  PID: 91515;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.886GB

Starting cleanup: 52 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:57:50
*  Total elapsed time (all runs):  1:10:48
*         Peak memory (this run):  3.5342 GB
*        Pipeline completed time: 2020-02-27 10:21:26
