### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name tutorial --genome hg38 --input /scratch/jps3dp/DATA/proseq/data/fastq/tutorial_r1.fq.gz --single-or-paired paired --input2 /scratch/jps3dp/DATA/proseq/data/fastq/tutorial_r2.fq.gz --protocol PROSEQ --prealignments human_rDNA -O /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline -M 8000`
*         Compute host:  udc-ba26-29c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/DATA/proseq
*            Outfolder:  /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/
*  Pipeline started at:   (11-27 13:55:04) elapsed: 5.0 _TIME_

### Version log:

*       Python version:  3.6.5
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.8.5
*        Pipeline hash:  1bce5897a1e47486f13b18f1f56420b279847613
*      Pipeline branch:  * dev
*        Pipeline date:  2019-11-27 11:21:24 -0500
*        Pipeline diff:  1 file changed, 1 insertion(+), 7 deletions(-)

### Arguments passed to pipeline:

*           `TSS_name`:  `None`
*            `adapter`:  `fastp`
*          `anno_name`:  `None`
*         `complexity`:  `False`
*        `config_file`:  `peppro.yaml`
*              `cores`:  `1`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `False`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/scratch/jps3dp/DATA/proseq/data/fastq/tutorial_r1.fq.gz']`
*             `input2`:  `['/scratch/jps3dp/DATA/proseq/data/fastq/tutorial_r2.fq.gz']`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `30`
*                `mem`:  `8000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline`
*         `paired_end`:  `True`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*           `protocol`:  `PROSEQ`
*            `recover`:  `False`
*        `sample_name`:  `tutorial`
*              `scale`:  `False`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `paired`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `0`
*          `verbosity`:  `None`

----------------------------------------

Local input file: /scratch/jps3dp/DATA/proseq/data/fastq/tutorial_r1.fq.gz
Local input file: /scratch/jps3dp/DATA/proseq/data/fastq/tutorial_r2.fq.gz

> `File_mb`	50.42	PEPPRO	_RES_

> `Read_type`	paired	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (11-27 13:55:05) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/tutorial_R1.fastq.gz`  

> `ln -sf /scratch/jps3dp/DATA/proseq/data/fastq/tutorial_r1.fq.gz /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/tutorial_R1.fastq.gz` (318626)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 318626;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/tutorial_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/tutorial_R2.fastq.gz`  

> `ln -sf /scratch/jps3dp/DATA/proseq/data/fastq/tutorial_r2.fq.gz /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/tutorial_R2.fastq.gz` (318628)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 318628;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/tutorial_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Missing stat 'Raw_reads'
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R1.fastq`,`/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2.fastq`  

> `gzip -f -d -c /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/tutorial_R1.fastq.gz > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R1.fastq` (318629)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.002GB.  
  PID: 318629;	Command: gzip;	Return code: 0;	Memory used: 0.002GB


> `gzip -f -d -c /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/tutorial_R2.fastq.gz > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2.fastq` (318634)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.002GB.  
  PID: 318634;	Command: gzip;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	2000000	PEPPRO	_RES_

> `Fastq_reads`	2000000	PEPPRO	_RES_
['/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/tutorial_R1.fastq.gz', '/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/tutorial_R2.fastq.gz']

### FASTQ processing:  (11-27 13:55:11) elapsed: 5.0 _TIME_

Missing stat 'Aligned_reads'
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq`  

> `( (fastp --overrepresentation_analysis --thread 1 --in1 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R1.fastq --adapter_sequence TGGAATTCTCGGGTGCCAAGG --length_required 18 --html /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastp/tutorial_R1_fastp_adapter.html --json /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastp/tutorial_R1_fastp_adapter.json --report_title 'tutorial' --stdout ) 2> /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastp/tutorial_R1_fastp_adapter.txt | seqtk trimfq -L 30 - | seqtk seq -r - > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq ) 2> /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastqc/tutorial_R1_rmAdapter.txt` (318658)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 0.123GB.  
  PID: 318658;	Command: ;	Return code: 0;	Memory used: 0.123GB

Evaluating read trimming

> `Trimmed_reads`	497796	PEPPRO	_RES_

> `Trim_loss_rate`	75.11	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastqc /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq` (318685)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of tutorial_R1_processed.fastq
Approx 5% complete for tutorial_R1_processed.fastq
Approx 10% complete for tutorial_R1_processed.fastq
Approx 15% complete for tutorial_R1_processed.fastq
Approx 20% complete for tutorial_R1_processed.fastq
Approx 25% complete for tutorial_R1_processed.fastq
Approx 30% complete for tutorial_R1_processed.fastq
Approx 35% complete for tutorial_R1_processed.fastq
Approx 40% complete for tutorial_R1_processed.fastq
Approx 45% complete for tutorial_R1_processed.fastq
Approx 50% complete for tutorial_R1_processed.fastq
Approx 55% complete for tutorial_R1_processed.fastq
Approx 60% complete for tutorial_R1_processed.fastq
Approx 65% complete for tutorial_R1_processed.fastq
Approx 70% complete for tutorial_R1_processed.fastq
Approx 75% complete for tutorial_R1_processed.fastq
Approx 80% complete for tutorial_R1_processed.fastq
Approx 85% complete for tutorial_R1_processed.fastq
Approx 90% complete for tutorial_R1_processed.fastq
Approx 95% complete for tutorial_R1_processed.fastq
Analysis complete for tutorial_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 0.164GB.  
  PID: 318685;	Command: fastqc;	Return code: 0;	Memory used: 0.164GB

> `FastQC report r1`	fastqc/tutorial_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq`  

> `( (fastp --overrepresentation_analysis --thread 1 --in1 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2.fastq --adapter_sequence GATCGTCGGACTGTAGAACTCTGAAC --length_required 18 --html /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastp/tutorial_R2_fastp_adapter.html --json /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastp/tutorial_R2_fastp_adapter.json --report_title 'tutorial' --stdout ) 2> /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastp/tutorial_R1_fastp_adapter.txt | seqtk trimfq -L 30 - | seqtk seq -r - > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq ) 2> /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastqc/tutorial_R1_rmAdapter.txt` (318711)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 0.429GB.  
  PID: 318711;	Command: ;	Return code: 0;	Memory used: 0.429GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq` (318724)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.429GB.  
  PID: 318724;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq.paired.fq`,`/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq.paired.fq`  

> `fastq_pair -t 1800000 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq` (318725)
<pre>
Left paired: 452944		Right paired: 452944
Left single: 44852		Right single: 20558
Writing the paired reads to /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq.paired.fq and /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq.single.fq and /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.429GB.  
  PID: 318725;	Command: fastq_pair;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/repaired.flag`  

> `mv /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq.paired.fq /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq` (318728)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.429GB.  
  PID: 318728;	Command: mv;	Return code: 0;	Memory used: 0.0GB


> `mv /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq.paired.fq /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq` (318730)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.429GB.  
  PID: 318730;	Command: mv;	Return code: 0;	Memory used: 0.0GB


> `touch repaired.flag` (318731)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.429GB.  
  PID: 318731;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (11-27 13:55:46) elapsed: 35.0 _TIME_

Skipping sample degradation plotting...
This requires using 'cutadapt' for adapter removal.

### Prealignments (11-27 13:55:46) elapsed: 0.0 _TIME_

Missing stat 'Aligned_reads'
Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (11-27 13:55:46) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/human_rDNA_bt2` (318732)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.429GB.  
  PID: 318732;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_summary.log`,`/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/human_rDNA_bt2 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R1.fq /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq` (318733)
<pre>
</pre>
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_summary.log`,`/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq.gz`  

> `(bowtie2 -p 1 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id tutorial -U /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq --un /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/human_rDNA_bt2 > /dev/null) 2>/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_summary.log` (318734)
<pre>
not gzipping output
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 0.429GB.  
  PID: 318734;	Command: bowtie2;	Return code: 0;	Memory used: 0.032GB


> `grep 'aligned exactly 1 time' /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_summary.log | awk '{print $1}'`

> `Aligned_reads_human_rDNA`	5860.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	1.18	PEPPRO	_RES_

### Map to genome (11-27 13:55:54) elapsed: 9.0 _TIME_

Missing stat 'Aligned_reads'
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam`  

> `bowtie2 -p 1 --very-sensitive -X 2000 --rg-id tutorial -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tmp7b3a497a -o /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam` (318749,318750,318751)
<pre>
2930 reads skipped
0 reads lost
450014 reads; of these:
  450014 (100.00%) were paired; of these:
    399197 (88.71%) aligned concordantly 0 times
    36510 (8.11%) aligned concordantly exactly 1 time
    14307 (3.18%) aligned concordantly >1 times
    ----
    399197 pairs aligned concordantly 0 times; of these:
      1862 (0.47%) aligned discordantly 1 time
    ----
    397335 pairs aligned 0 times concordantly or discordantly; of these:
      794670 mates make up the pairs; of these:
        468224 (58.92%) aligned 0 times
        239606 (30.15%) aligned exactly 1 time
        86840 (10.93%) aligned >1 times
47.98% overall alignment rate
</pre>
Command completed. Elapsed time: 0:05:51. Running peak memory: 3.465GB.  
  PID: 318749;	Command: bowtie2;	Return code: 0;	Memory used: 3.465GB  
  PID: 318751;	Command: samtools;	Return code: 0;	Memory used: 0.201GB  
  PID: 318750;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


> `samtools view -q 10 -b -@ 1 -U /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_fail_qc.bam /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam` (319341)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.465GB.  
  PID: 319341;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	431804	PEPPRO	_RES_

> `QC_filtered_reads`	274898	PEPPRO	_RES_

> `Aligned_reads`	156906.0	PEPPRO	_RES_

> `Alignment_rate`	31.52	PEPPRO	_RES_

> `Total_efficiency`	7.85	PEPPRO	_RES_

> `Read_depth`	1.17	PEPPRO	_RES_

### Compress all unmapped read files (11-27 14:02:12) elapsed: 378.0 _TIME_

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R1.fq.gz`  

> `gzip -f /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R1.fq` (319372)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.465GB.  
  PID: 319372;	Command: gzip;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq.gz`  

> `gzip -f /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq` (319378)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.465GB.  
  PID: 319378;	Command: gzip;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam.bai`  

> `samtools index /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam` (319386)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319386;	Command: samtools;	Return code: 0;	Memory used: 0.001GB


> `samtools idxstats /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	5357	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_noMT.bam`  

> `samtools index /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam` (319392)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319392;	Command: samtools;	Return code: 0;	Memory used: 0.0GB


> `samtools idxstats /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 1 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_noMT.bam` (319394,319395,319396,319397)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.465GB.  
  PID: 319394;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 319396;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 319395;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 319397;	Command: xargs;	Return code: 0;	Memory used: 0.011GB


> `mv /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_noMT.bam /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam` (319402)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319402;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam` (319403)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319403;	Command: samtools;	Return code: 0;	Memory used: 0.0GB


### Split BAM file (11-27 14:02:25) elapsed: 13.0 _TIME_

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam`,`/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam | samtools sort - -@ 1 > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam` (319406,319407)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319406;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 319407;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools view -b -f 128 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam | samtools sort - -@ 1 > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE2.bam` (319411,319412)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.465GB.  
  PID: 319411;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 319412;	Command: samtools;	Return code: 0;	Memory used: 0.068GB

Missing stat 'Maximum_read_length'

> `Maximum_read_length`	30	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (11-27 14:02:30) elapsed: 5.0 _TIME_

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam` (319418)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319418;	Command: samtools;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam -c 1 -o /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_bamQC.tsv` (319419)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tmp_tutorial_PE1_hfixb9vw'
Processing with 1 cores...
Discarding 158 chunk(s) of reads: ['chrM', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2']
Keeping 37 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270708v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_KI270725v1_random', 'chr17_GL000205v2_random', 'chr22_KI270733v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_KI270750v1', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.465GB.  
  PID: 319419;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.034GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	15975.0	PEPPRO	_RES_

> `PBC2`	15975.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_unmap.bam`  

> `samtools view -b -@ 1 -f 12  /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_unmap.bam` (319425)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.465GB.  
  PID: 319425;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -c -f 4 -@ 1 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam`

> `Unmapped_reads`	468224	PEPPRO	_RES_

### Split BAM by strand (11-27 14:02:35) elapsed: 5.0 _TIME_

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam`,`/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam` (319430)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319430;	Command: samtools;	Return code: 0;	Memory used: 0.002GB


> `samtools view -bh -f 16 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam` (319432)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319432;	Command: samtools;	Return code: 0;	Memory used: 0.002GB


### Calculate TSS enrichment (11-27 14:02:36) elapsed: 1.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (319434)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319434;	Command: sed;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/plus_TSS.tsv -p ends -c 1 -z -v -s 6 -o /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_TssEnrichment.txt` (319436)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.465GB.  
  PID: 319436;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.078GB


> `TSS_Plus_Score`	33.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/minus_TSS.tsv -p ends -c 1 -z -v -s 6 -o /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_TssEnrichment.txt` (319444)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.465GB.  
  PID: 319444;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.079GB


> `TSS_Minus_Score`	4.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_TssEnrichment.txt /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_TssEnrichment.txt` (319451)
<pre>

Generating TSS plot with /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_TssEnrichment.txt and /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.465GB.  
  PID: 319451;	Command: Rscript;	Return code: 0;	Memory used: 0.136GB

> `TSS enrichment`	QC_hg38/tutorial_TSSenrichment.pdf	TSS enrichment	QC_hg38/tutorial_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt` (319486,319488,319489,319490)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319486;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 319489;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 319488;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 319490;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt` (319492)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319492;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (11-27 14:02:51) elapsed: 16.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_tss.bed` (319494,319495)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.465GB.  
  PID: 319494;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 319495;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_gene_body.bed` (319499,319500)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319499;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 319500;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_TSS_density.bed` (319503,319504,319505,319506)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319503;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 319505;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 319504;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 319506;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_gene_body_density.bed` (319508,319509,319510)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319509;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 319508;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB  
  PID: 319510;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_TSS_density.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.bed` (319513,319514,319515)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319513;	Command: join;	Return code: 0;	Memory used: 0.0GB  
  PID: 319514;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 319515;	Command: env;	Return code: 0;	Memory used: 0.003GB


> `sort -k5,5n /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	228.97	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi -i /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.bed` (319521)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.465GB.  
  PID: 319521;	Command: Rscript;	Return code: 0;	Memory used: 0.201GB

> `Pause index`	QC_hg38/tutorial_pause_index.pdf	Pause index	QC_hg38/tutorial_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.bed.gz`  

> `gzip -f -f /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.bed` (319553)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319553;	Command: gzip;	Return code: 0;	Memory used: 0.0GB


### Calculate FRiP (11-27 14:03:04) elapsed: 12.0 _TIME_


> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam`
156906.0 11064

> `Plus FRiP`	0.07	PEPPRO	_RES_

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam`
156906.0 10420

> `Minus FRiP`	0.07	PEPPRO	_RES_

### Plot fragment distribution (11-27 14:03:04) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLenDistribution.pdf`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/fragment_length_dist.pl /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLen.txt` (319566)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319566;	Command: perl;	Return code: 0;	Memory used: 0.0GB


> `sort -n  /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLen.txt | uniq -c  > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_fragCount.txt` (319568,319569)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319568;	Command: sort;	Return code: 0;	Memory used: 0.0GB  
  PID: 319569;	Command: uniq;	Return code: 0;	Memory used: 0.002GB


> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frag -l /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLen.txt -c /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_fragCount.txt -p /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLenDistribution.pdf -t /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLenDistribution.txt` (319571)
<pre>
Fragment distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.465GB.  
  PID: 319571;	Command: Rscript;	Return code: 0;	Memory used: 0.234GB

> `Fragment distribution`	QC_hg38/tutorial_fragLenDistribution.pdf	Fragment distribution	QC_hg38/tutorial_fragLenDistribution.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/hg38_annotations.bed.gz` (319602)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319602;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `gzip -f -d -c /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/hg38_annotations.bed.gz > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/hg38_annotations.bed` (319603)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319603;	Command: gzip;	Return code: 0;	Memory used: 0.002GB


### Calculate fraction of reads in features (FRiF) (11-27 14:03:10) elapsed: 6.0 _TIME_


> `cut -f 4 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/"$4}' /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/raw/hg38_annotations.bed` (319609)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319609;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/3' UTR" "/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/3_UTR"` (319612)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319612;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/3_UTR_sort.bed` (319613,319614,319615,319616)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319613;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 319614;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 319616;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 319615;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted -counts -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_plus_coverage.bed` (319618)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319618;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted -counts -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_minus_coverage.bed` (319621)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319621;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB

Target exists: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/5' UTR" "/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/5_UTR"` (319623)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319623;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/5_UTR_sort.bed` (319624,319625,319626,319627)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319624;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 319625;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 319627;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 319626;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted -counts -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_plus_coverage.bed` (319630)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319630;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted -counts -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_minus_coverage.bed` (319632)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319632;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB

Target exists: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Enhancer_sort.bed` (319634,319635,319636,319637)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319634;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 319636;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 319635;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 319637;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted -counts -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_plus_coverage.bed` (319640)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319640;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted -counts -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_minus_coverage.bed` (319642)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319642;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB

Target exists: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Exon_sort.bed` (319644,319645,319646,319647)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.465GB.  
  PID: 319644;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 319645;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 319647;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 319646;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted -counts -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_plus_coverage.bed` (319661)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319661;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted -counts -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_minus_coverage.bed` (319664)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319664;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB

Target exists: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Intron_sort.bed` (319667,319668,319669,319670)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.465GB.  
  PID: 319667;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 319670;	Command: bedtools;	Return code: 0;	Memory used: 0.073GB  
  PID: 319668;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 319669;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted -counts -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_plus_coverage.bed` (319704)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319704;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted -counts -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_minus_coverage.bed` (319706)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319706;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB

Target exists: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter_sort.bed` (319709,319710,319711,319712)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319709;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 319710;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 319712;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 319711;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted -counts -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_plus_coverage.bed` (319714)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319714;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted -counts -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_minus_coverage.bed` (319716)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319716;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB

Target exists: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region"` (319718)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319718;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region_sort.bed` (319719,319721,319722,319723)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.465GB.  
  PID: 319719;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 319722;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 319721;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 319723;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted -counts -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_plus_coverage.bed` (319726)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319726;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted -counts -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_minus_coverage.bed` (319728)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319728;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB


### Plot FRiF (11-27 14:03:29) elapsed: 19.0 _TIME_


> `samtools view -@ 1 -q 10 -c -F4 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam`
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -n tutorial -r 16265 -o /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_frif.pdf --bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_plus_coverage.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_plus_coverage.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_plus_coverage.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_plus_coverage.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_plus_coverage.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_plus_coverage.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_plus_coverage.bed` (319731)
<pre>
Cumulative FRiF plot completed!

</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.465GB.  
  PID: 319731;	Command: Rscript;	Return code: 0;	Memory used: 0.44GB

> `Plus FRiF`	QC_hg38/tutorial_plus_frif.pdf	Plus FRiF	QC_hg38/tutorial_plus_frif.png	PEPPRO	_OBJ_

> `samtools view -@ 1 -q 10 -c -F4 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam`
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -n tutorial -r 15685 -o /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_frif.pdf --bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_minus_coverage.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_minus_coverage.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_minus_coverage.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_minus_coverage.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_minus_coverage.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_minus_coverage.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_minus_coverage.bed` (319781)
<pre>
Cumulative FRiF plot completed!

</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.465GB.  
  PID: 319781;	Command: Rscript;	Return code: 0;	Memory used: 0.446GB

> `Minus FRiF`	QC_hg38/tutorial_minus_frif.pdf	Minus FRiF	QC_hg38/tutorial_minus_frif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (11-27 14:04:23) elapsed: 54.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/hg38_exons_sort.bed` (319833,319834)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.465GB.  
  PID: 319834;	Command: bedtools;	Return code: 0;	Memory used: 0.08GB  
  PID: 319833;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/hg38_introns_sort.bed` (319841,319842,319843)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.465GB.  
  PID: 319842;	Command: bedtools;	Return code: 0;	Memory used: 0.084GB  
  PID: 319841;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 319843;	Command: bedtools;	Return code: 0;	Memory used: 0.092GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_coverage.bed`,`/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_coverage.bed` (319851)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319851;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam -g /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_coverage.bed` (319853)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319853;	Command: bedtools;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/0.156906)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_rpkm.bed` (319862,319863,319864)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319862;	Command: awk;	Return code: 0;	Memory used: 0.009GB  
  PID: 319864;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 319863;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/0.156906)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_rpkm.bed` (319866,319867,319868)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.465GB.  
  PID: 319866;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 319868;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 319867;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_rpkm.bed /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_exon_intron_ratios.bed` (319871,319872,319873)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319871;	Command: join;	Return code: 0;	Memory used: 0.0GB  
  PID: 319872;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 319873;	Command: sort;	Return code: 0;	Memory used: 0.002GB


> `awk '{print $5}' /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	4.26	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_exon_intron_ratios.bed --raw` (319879)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.465GB.  
  PID: 319879;	Command: Rscript;	Return code: 0;	Memory used: 0.215GB

> `mRNA contamination`	QC_hg38/tutorial_mRNA_contamination.pdf	mRNA contamination	QC_hg38/tutorial_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_exon_intron_ratios.bed.gz`  

> `gzip -f -f /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/QC_hg38/tutorial_exon_intron_ratios.bed` (319922)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319922;	Command: gzip;	Return code: 0;	Memory used: 0.001GB


### Produce bigWig files (11-27 14:04:44) elapsed: 21.0 _TIME_

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/signal_hg38/tutorial_plus_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam` (319924)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 319924;	Command: samtools;	Return code: 0;	Memory used: 0.0GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/signal_hg38/tutorial_plus_body_0-mer.bw -p 1 --variable-step --tail-edge` (319925)
<pre>
Registering input file: '/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam'
Temporary files will be stored in: 'tmp_tutorial_plus_cuttrace_oxinje4r'
Processing with 1 cores...
Discarding 165 chunk(s) of reads: ['chrM', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Keeping 30 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr14_KI270722v1_random', 'chr14_KI270725v1_random', 'chr22_KI270733v1_random', 'chr22_KI270738v1_random', 'chrUn_GL000195v1', 'chrUn_GL000219v1']
Reduce step (merge files)...
Merging 30 files into output file: '/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/signal_hg38/tutorial_plus_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.465GB.  
  PID: 319925;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.032GB

Target to produce: `/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/signal_hg38/tutorial_minus_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam` (320065)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.465GB.  
  PID: 320065;	Command: samtools;	Return code: 0;	Memory used: 0.001GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/signal_hg38/tutorial_minus_body_0-mer.bw -p 1 --variable-step --tail-edge` (320066)
<pre>
Registering input file: '/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam'
Temporary files will be stored in: 'tmp_tutorial_minus_cuttrace_785gqf3g'
Processing with 1 cores...
Discarding 162 chunk(s) of reads: ['chrM', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2']
Keeping 33 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270708v1_random', 'chr14_GL000225v1_random', 'chr17_GL000205v2_random', 'chrUn_KI270442v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_KI270750v1', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 33 files into output file: '/project/shefflab/processed/proseq/peppro_tutorial/pe/11-27-19/results_pipeline/tutorial/signal_hg38/tutorial_minus_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.465GB.  
  PID: 320066;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.032GB

Starting cleanup: 55 files; 3 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:09:52
*  Total elapsed time (all runs):  0:21:54
*         Peak memory (this run):  3.4646 GB
*        Pipeline completed time: 2019-11-27 14:04:51
