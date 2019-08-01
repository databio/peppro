### [Pipeline run code and environment:]

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name tutorial --genome hg38 --input /scratch/jps3dp/DATA/proseq/tutorial_r1.fq.gz --single-or-paired paired --input2 /scratch/jps3dp/DATA/proseq/tutorial_r2.fq.gz --prealignments human_rDNA rCRSd -O /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline -M 8000`
*         Compute host:  udc-ba26-32c1
*          Working dir:  /sfs/lustre/scratch/jps3dp/DATA/proseq
*            Outfolder:  /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/
*  Pipeline started at:   (07-31 20:43:17) elapsed: 1.0 _TIME_

### [Version log:]

*       Python version:  3.6.5
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.11.3
*         Pipeline dir:  `/sfs/lustre/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.8.0
*        Pipeline hash:  a9b81015dd8bff626b2aa7be9e9038fa5ff2937a
*      Pipeline branch:  * dev
*        Pipeline date:  2019-07-31 14:14:51 -0400
*        Pipeline diff:  5 files changed, 44 insertions(+), 44 deletions(-)

### [Arguments passed to pipeline:]

*            `recover`:  `False`
*          `new_start`:  `False`
*              `dirty`:  `False`
*       `force_follow`:  `False`
*        `config_file`:  `peppro.yaml`
*      `output_parent`:  `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline`
*                `mem`:  `8000`
*              `cores`:  `1`
*        `sample_name`:  `tutorial`
*              `input`:  `['/scratch/jps3dp/DATA/proseq/tutorial_r1.fq.gz']`
*             `input2`:  `['/scratch/jps3dp/DATA/proseq/tutorial_r2.fq.gz']`
*    `genome_assembly`:  `hg38`
*   `single_or_paired`:  `paired`
*              `runon`:  `pro`
*            `adapter`:  `fastp`
*              `dedup`:  `seqkit`
*            `trimmer`:  `seqtk`
*                `umi`:  `False`
*            `umi_len`:  `8`
*            `max_len`:  `30`
*                `sob`:  `False`
*              `scale`:  `False`
*              `parts`:  `4`
*      `prealignments`:  `['human_rDNA', 'rCRSd']`
*           `TSS_name`:  `None`
*        `ensembl_tss`:  `None`
*  `ensembl_gene_body`:  `None`
*           `pre_name`:  `None`
*          `anno_name`:  `None`
*          `exon_name`:  `None`
*        `intron_name`:  `None`
*           `coverage`:  `False`
*               `keep`:  `False`
*            `no_fifo`:  `False`
*         `complexity`:  `True`
*         `paired_end`:  `True`

----------------------------------------


Changed status from initializing to running.

Loading config file: /scratch/jps3dp/tools/databio/peppro/pipelines/peppro.yaml

Local input file: /scratch/jps3dp/DATA/proseq/tutorial_r1.fq.gz
Local input file: /scratch/jps3dp/DATA/proseq/tutorial_r2.fq.gz

> `File_mb`	50.42	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `Read_type`	paired	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `Genome`	hg38	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

### Merge/link and fastq conversion:  (07-31 20:43:18) elapsed: 1.0 _TIME_

Number of input file sets: 2
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/tutorial_R1.fastq.gz']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/tutorial_R1.fastq.gz`


> `ln -sf /scratch/jps3dp/DATA/proseq/tutorial_r1.fq.gz /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/tutorial_R1.fastq.gz` (418166)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.
  PID: 418166;	Command: ln;	Return code: 0;	Memory used: 0.0GB
Local input file: '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/tutorial_R1.fastq.gz'
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/tutorial_R2.fastq.gz']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/tutorial_R2.fastq.gz`


> `ln -sf /scratch/jps3dp/DATA/proseq/tutorial_r2.fq.gz /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/tutorial_R2.fastq.gz` (418167)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.
  PID: 418167;	Command: ln;	Return code: 0;	Memory used: 0.0GB
Local input file: '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/tutorial_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1.fastq', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2.fastq']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1.fastq`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2.fastq`


> `gzip -f -d -c /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/tutorial_R1.fastq.gz > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1.fastq` (418169)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.001GB.
  PID: 418169;	Command: gzip;	Return code: 0;	Memory used: 0.001GB

> `gzip -f -d -c /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/tutorial_R2.fastq.gz > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2.fastq` (418173)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.001GB.
  PID: 418173;	Command: gzip;	Return code: 0;	Memory used: 0.001GB
Follow:

> `Raw_reads`	2000000	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `Fastq_reads`	2000000	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv
['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/tutorial_R1.fastq.gz', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/tutorial_R2.fastq.gz']

### FASTQ processing:  (07-31 20:43:22) elapsed: 4.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq`


> `(fastp --overrepresentation_analysis --thread 1 --in1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1.fastq --adapter_sequence TGGAATTCTCGGGTGCCAAGG --length_required 26 --html /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastqc/tutorial_R1_rmAdapter.html --json /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastqc/tutorial_R1_rmAdapter.json --report_title 'tutorial' -o /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_noadap.fastq ) 2> /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastqc/tutorial_R1_rmAdapter.txt` (418211)

<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 0.535GB.
  PID: 418211;	Command: fastp;	Return code: 0;	Memory used: 0.535GB

> `seqkit rmdup --threads 1 --by-seq --ignore-case -o /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_dedup.fastq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_noadap.fastq` (418223)

<pre>
[INFO][0m 3022 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.535GB.
  PID: 418223;	Command: seqkit;	Return code: 0;	Memory used: 0.0GB

> `seqtk trimfq -b 8 -L 30 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_noadap.fastq | seqtk seq -r - > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq` (418234,418241)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.535GB.
  PID: 418234;	Command: seqtk;	Return code: 0;	Memory used: 0.0GB
  PID: 418241;	Command: seqtk;	Return code: 0;	Memory used: 0.0GB
Follow:
> `FastP_report`	fastqc/tutorial_R1_rmAdapter.html	FastP_report	None	PEPPRO	_OBJ_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/objects.tsv

> `grep 'reads with adapter trimmed:' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastqc/tutorial_R1_rmAdapter.txt | head -n 1 | awk '{print $NF}'`


> `Reads_with_adapter`	598996.0	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `grep 'total bases:' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastqc/tutorial_R1_rmAdapter.txt | head -n 1 | awk '{print $NF}'`


> `grep 'bases trimmed due to adapters:' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastqc/tutorial_R1_rmAdapter.txt | awk '{print $NF}'`


> `Pct_adapter_contamination`	0.44	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `grep 'reads failed due to too short:' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastqc/tutorial_R1_rmAdapter.txt | head -n 1 | awk '{print $NF}'`


> `Reads_too_short`	522549.0	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `Duplicate_reads`	3022.0	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq`


> `seqtk trimfq -b 8 -L 30 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_dedup.fastq | seqtk seq -r - > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq` (418276,418277)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.535GB.
  PID: 418276;	Command: seqtk;	Return code: 0;	Memory used: 0.0GB
  PID: 418277;	Command: seqtk;	Return code: 0;	Memory used: 0.0GB
Follow:
Evaluating read trimming

> `Trimmed_reads`	461368	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `Trim_loss_rate`	76.93	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq`


> `(fastp --overrepresentation_analysis --thread 1 --in1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2.fastq --adapter_sequence GATCGTCGGACTGTAGAACTCTGAAC --length_required 26 --html /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastqc/tutorial_R2_rmAdapter.html --json /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastqc/tutorial_R2_rmAdapter.json --report_title 'tutorial' --stdout ) 2> /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastqc/tutorial_R1_rmAdapter.txt | seqtk trimfq -e 8 - | seqtk seq -r - > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq` (418284,418285,418286)

<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.535GB.
  PID: 418285;	Command: seqtk;	Return code: 0;	Memory used: 0.0GB
  PID: 418284;	Command: fastp;	Return code: 0;	Memory used: 0.091GB
  PID: 418286;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq`


> `cp /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq` (418296)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.535GB.
  PID: 418296;	Command: cp;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq.paired.fq', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq.paired.fq']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq.paired.fq`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq.paired.fq`


> `fastq_pair -t 1800000 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq` (418298)

<pre>
Left paired: 418940		Right paired: 418940
Left single: 42428		Right single: 20579
Writing the paired reads to /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq.paired.fq and /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq.paired.fq.
Writing the single reads to /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq.single.fq and /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.535GB.
  PID: 418298;	Command: fastq_pair;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/repaired.flag']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/repaired.flag`


> `mv /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq.paired.fq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq` (418301)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.535GB.
  PID: 418301;	Command: mv;	Return code: 0;	Memory used: 0.0GB

> `mv /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq.paired.fq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq` (418302)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.535GB.
  PID: 418302;	Command: mv;	Return code: 0;	Memory used: 0.0GB

> `touch repaired.flag` (418303)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.535GB.
  PID: 418303;	Command: touch;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq.paired.fq', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq.paired.fq']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq.paired.fq`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq.paired.fq`


> `fastq_pair -t 1800000 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq` (418304)

<pre>
Left paired: 420772		Right paired: 420772
Left single: 43618		Right single: 18747
Writing the paired reads to /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq.paired.fq and /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq.single.fq and /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.535GB.
  PID: 418304;	Command: fastq_pair;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/dups_repaired.flag']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/dups_repaired.flag`


> `mv /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq.paired.fq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq` (418306)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.535GB.
  PID: 418306;	Command: mv;	Return code: 0;	Memory used: 0.0GB

> `mv /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq.paired.fq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq` (418307)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.535GB.
  PID: 418307;	Command: mv;	Return code: 0;	Memory used: 0.0GB

> `touch dups_repaired.flag` (418308)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.535GB.
  PID: 418308;	Command: touch;	Return code: 0;	Memory used: 0.0GB

### Prealignments (07-31 20:43:43) elapsed: 21.0 _TIME_

Prealignment assemblies: ['human_rDNA', 'rCRSd']

### Map to human_rDNA (07-31 20:43:43) elapsed: 0.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/human_rDNA_bt2']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/human_rDNA_bt2`


> `mkfifo /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/human_rDNA_bt2` (418309)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.535GB.
  PID: 418309;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_summary.log', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq.gz']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_summary.log`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq.gz`


> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/human_rDNA_bt2 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R1.fq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq` (418311)

<pre>
</pre>
Not waiting for subprocesses: [418311]
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_summary.log', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq.gz']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_summary.log`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq.gz`


> `(bowtie2 -p 1 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /nv/t1/genomes/human_rDNA/bowtie2_index/human_rDNA --rg-id tutorial -U /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq --un /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/human_rDNA_bt2 > /dev/null) 2>/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_summary.log` (418312)

<pre>
not gzipping output
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.535GB.
  PID: 418312;	Command: bowtie2;	Return code: 0;	Memory used: 0.022GB

> `grep 'aligned exactly 1 time' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_summary.log | awk '{print $1}'`


> `Aligned_reads_human_rDNA`	68452.0	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `Alignment_rate_human_rDNA`	14.84	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

### Map to human_rDNA (07-31 20:43:51) elapsed: 8.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/human_rDNA_dups_bt2']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/human_rDNA_dups_bt2`


> `mkfifo /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/human_rDNA_dups_bt2` (418446)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.535GB.
  PID: 418446;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_dups_summary.log', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_dups_R2.fq.gz']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_dups_summary.log`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_dups_R2.fq.gz`


> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/human_rDNA_dups_bt2 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_dups_R1.fq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_dups_R2.fq` (418447)

<pre>
</pre>
Not waiting for subprocesses: [418447]
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_dups_summary.log', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_dups_R2.fq.gz']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_dups_summary.log`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_dups_R2.fq.gz`


> `(bowtie2 -p 1 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /nv/t1/genomes/human_rDNA/bowtie2_index/human_rDNA --rg-id tutorial -U /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq --un /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/human_rDNA_dups_bt2 > /dev/null) 2>/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_bt_aln_dups_summary.log` (418448)

<pre>
not gzipping output
34226 reads skipped
0 reads lost
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 0.535GB.
  PID: 418448;	Command: bowtie2;	Return code: 0;	Memory used: 0.022GB

### Map to rCRSd (07-31 20:43:59) elapsed: 8.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/rCRSd_bt2']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/rCRSd_bt2`


> `mkfifo /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/rCRSd_bt2` (418461)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.535GB.
  PID: 418461;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_bt_aln_summary.log', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_R2.fq.gz']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_bt_aln_summary.log`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_R2.fq.gz`


> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/rCRSd_bt2 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R1.fq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_R1.fq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_R2.fq` (418462)

<pre>
</pre>
Not waiting for subprocesses: [418462]
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_bt_aln_summary.log', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_R2.fq.gz']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_bt_aln_summary.log`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_R2.fq.gz`


> `(bowtie2 -p 1 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /nv/t1/genomes/rCRSd/bowtie2_index/rCRSd --rg-id tutorial -U /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R1.fq --un /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/rCRSd_bt2 > /dev/null) 2>/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_bt_aln_summary.log` (418463)

<pre>
not gzipping output
34533 reads skipped
0 reads lost
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.535GB.
  PID: 418463;	Command: bowtie2;	Return code: 0;	Memory used: 0.024GB

> `grep 'aligned exactly 1 time' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_bt_aln_summary.log | awk '{print $1}'`


> `Aligned_reads_rCRSd`	12998.0	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `Alignment_rate_rCRSd`	2.82	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

### Map to rCRSd (07-31 20:44:08) elapsed: 8.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/rCRSd_dups_bt2']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/rCRSd_dups_bt2`


> `mkfifo /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/rCRSd_dups_bt2` (418497)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.535GB.
  PID: 418497;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_bt_aln_dups_summary.log', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_dups_R2.fq.gz']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_bt_aln_dups_summary.log`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_dups_R2.fq.gz`


> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/rCRSd_dups_bt2 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_dups_R1.fq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_dups_R2.fq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_dups_R1.fq /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_dups_R2.fq` (418498)

<pre>
</pre>
Not waiting for subprocesses: [418498]
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_bt_aln_dups_summary.log', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_dups_R2.fq.gz']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_bt_aln_dups_summary.log`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_dups_R2.fq.gz`


> `(bowtie2 -p 1 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /nv/t1/genomes/rCRSd/bowtie2_index/rCRSd --rg-id tutorial -U /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_dups_R1.fq --un /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/rCRSd_dups_bt2 > /dev/null) 2>/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_bt_aln_dups_summary.log` (418499)

<pre>
not gzipping output
6499 reads skipped
0 reads lost
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.535GB.
  PID: 418499;	Command: bowtie2;	Return code: 0;	Memory used: 0.024GB

### Map to genome (07-31 20:44:15) elapsed: 7.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam`

6503 reads skipped
0 reads lost

> `bowtie2 -p 1 --very-sensitive -X 2000 --rg-id tutorial -x /nv/t1/genomes/hg38/bowtie2_index/hg38 --rf -1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_R1.fq -2 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tmplh2_a83a -o /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam` (418519,418530,418532)

<pre>
378215 reads; of these:
  378215 (100.00%) were paired; of these:
    41435 (10.96%) aligned concordantly 0 times
    279989 (74.03%) aligned concordantly exactly 1 time
    56791 (15.02%) aligned concordantly >1 times
    ----
    41435 pairs aligned concordantly 0 times; of these:
      8119 (19.59%) aligned discordantly 1 time
    ----
    33316 pairs aligned 0 times concordantly or discordantly; of these:
      66632 mates make up the pairs; of these:
        44906 (67.39%) aligned 0 times
        13697 (20.56%) aligned exactly 1 time
        8029 (12.05%) aligned >1 times
94.06% overall alignment rate
</pre>
Command completed. Elapsed time: 0:04:30. Running peak memory: 3.45GB.
  PID: 418532;	Command: samtools;	Return code: 0;	Memory used: 0.0GB
  PID: 418519;	Command: bowtie2;	Return code: 0;	Memory used: 3.45GB
  PID: 418530;	Command: samtools;	Return code: 0;	Memory used: 0.0GB

> `samtools view -q 10 -b -@ 1 -U /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_fail_qc.bam /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam` (419432)

<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.45GB.
  PID: 419432;	Command: samtools;	Return code: 0;	Memory used: 0.01GB
Follow:

> `samtools depth -b /nv/t1/genomes/hg38/refgene_anno/hg38_pre-mRNA.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`


> `Mapped_reads`	711524	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `QC_filtered_reads`	384789	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `Aligned_reads`	326735.0	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `Alignment_rate`	70.82	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `Total_efficiency`	16.34	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `Read_depth`	1.35	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort_dups.bam']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort_dups.bam`


> `bowtie2 -p 1 --very-sensitive -X 2000 --rg-id tutorial -x /nv/t1/genomes/hg38/bowtie2_index/hg38 --rf -1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_dups_R1.fq -2 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tmplh2_a83a -o /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp_dups.bam` (419470,419471,419472)

<pre>
379736 reads; of these:
  379736 (100.00%) were paired; of these:
    41728 (10.99%) aligned concordantly 0 times
    280272 (73.81%) aligned concordantly exactly 1 time
    57736 (15.20%) aligned concordantly >1 times
    ----
    41728 pairs aligned concordantly 0 times; of these:
      8143 (19.51%) aligned discordantly 1 time
    ----
    33585 pairs aligned 0 times concordantly or discordantly; of these:
      67170 mates make up the pairs; of these:
        45195 (67.28%) aligned 0 times
        13717 (20.42%) aligned exactly 1 time
        8258 (12.29%) aligned >1 times
94.05% overall alignment rate
</pre>
Command completed. Elapsed time: 0:04:18. Running peak memory: 3.45GB.
  PID: 419471;	Command: samtools;	Return code: 0;	Memory used: 0.0GB
  PID: 419470;	Command: bowtie2;	Return code: 0;	Memory used: 3.45GB
  PID: 419472;	Command: samtools;	Return code: 0;	Memory used: 0.0GB

> `samtools view -q 10 -b -@ 1 -U /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_fail_qc_dups.bam /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp_dups.bam > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort_dups.bam` (419887)

<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.45GB.
  PID: 419887;	Command: samtools;	Return code: 0;	Memory used: 0.01GB

### Compress all unmapped read files (07-31 20:53:45) elapsed: 570.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R1.fq.gz']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R1.fq.gz`


> `gzip -f /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R1.fq` (419893)

<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.45GB.
  PID: 419893;	Command: gzip;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq.gz']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq.gz`


> `gzip -f /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_R2.fq` (419896)

<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.45GB.
  PID: 419896;	Command: gzip;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_R1.fq.gz']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_R1.fq.gz`


> `gzip -f /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_R1.fq` (419900)

<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.45GB.
  PID: 419900;	Command: gzip;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_R2.fq.gz']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_R2.fq.gz`


> `gzip -f /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_R2.fq` (419905)

<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.45GB.
  PID: 419905;	Command: gzip;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam.bai']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam.bai`


> `samtools index /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam` (419909)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 419909;	Command: samtools;	Return code: 0;	Memory used: 0.0GB

> `samtools idxstats /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`


> `Mitochondrial_reads`	249	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_noMT.bam']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_noMT.bam`


> `samtools index /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam` (419914)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 419914;	Command: samtools;	Return code: 0;	Memory used: 0.0GB

> `samtools idxstats /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_noMT.bam` (419916,419917,419918,419919)

<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.45GB.
  PID: 419916;	Command: samtools;	Return code: 0;	Memory used: 0.0GB
  PID: 419918;	Command: grep;	Return code: 0;	Memory used: 0.0GB
  PID: 419917;	Command: cut;	Return code: 0;	Memory used: 0.0GB
  PID: 419919;	Command: xargs;	Return code: 0;	Memory used: 0.0GB

> `mv /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_noMT.bam /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam` (419925)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 419925;	Command: mv;	Return code: 0;	Memory used: 0.0GB

> `samtools index /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam` (419926)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 419926;	Command: samtools;	Return code: 0;	Memory used: 0.0GB

> `samtools stats /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`


### Split BAM file (07-31 20:54:07) elapsed: 22.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE2.bam']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE2.bam`


> `samtools view -b -f 64 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam | samtools sort - -@ 1 > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam` (419935,419936)

<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.45GB.
  PID: 419935;	Command: samtools;	Return code: 0;	Memory used: 0.0GB
  PID: 419936;	Command: samtools;	Return code: 0;	Memory used: 0.087GB

> `samtools view -b -f 128 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam | samtools sort - -@ 1 > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE2.bam` (419946,419947)

<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.45GB.
  PID: 419946;	Command: samtools;	Return code: 0;	Memory used: 0.0GB
  PID: 419947;	Command: samtools;	Return code: 0;	Memory used: 0.082GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp_dups.bam.bai']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp_dups.bam.bai`


> `samtools index /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp_dups.bam` (419953)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 419953;	Command: samtools;	Return code: 0;	Memory used: 0.0GB

> `samtools idxstats /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_noMT_dups.bam']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_noMT_dups.bam`


> `samtools index /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort_dups.bam` (419959)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 419959;	Command: samtools;	Return code: 0;	Memory used: 0.0GB

> `samtools idxstats /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort_dups.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort_dups.bam > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_noMT_dups.bam` (419961,419962,419963,419964)

<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.45GB.
  PID: 419963;	Command: grep;	Return code: 0;	Memory used: 0.0GB
  PID: 419961;	Command: samtools;	Return code: 0;	Memory used: 0.0GB
  PID: 419962;	Command: cut;	Return code: 0;	Memory used: 0.0GB
  PID: 419964;	Command: xargs;	Return code: 0;	Memory used: 0.013GB

> `mv /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_noMT_dups.bam /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort_dups.bam` (419970)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 419970;	Command: mv;	Return code: 0;	Memory used: 0.0GB

> `samtools index /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort_dups.bam` (419971)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 419971;	Command: samtools;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_dups_PE1.bam', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_dups_PE2.bam']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_dups_PE1.bam`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_dups_PE2.bam`


> `samtools view -b -f 64 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort_dups.bam | samtools sort - -@ 1 > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_dups_PE1.bam` (419973,419974)

<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.45GB.
  PID: 419973;	Command: samtools;	Return code: 0;	Memory used: 0.0GB
  PID: 419974;	Command: samtools;	Return code: 0;	Memory used: 0.083GB

> `samtools view -b -f 128 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort_dups.bam | samtools sort - -@ 1 > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_dups_PE2.bam` (419981,419982)

<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.45GB.
  PID: 419981;	Command: samtools;	Return code: 0;	Memory used: 0.0GB
  PID: 419982;	Command: samtools;	Return code: 0;	Memory used: 0.08GB

### Calculate library complexity (07-31 20:54:35) elapsed: 28.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_out.txt']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_out.txt`


> `preseq c_curve -v -o /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_out.txt -B /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_dups_PE1.bam` (419989)

<pre>
BAM_INPUT
TOTAL READS     = 330356
COUNTS_SUM      = 330356
DISTINCT READS  = 324913
DISTINCT COUNTS = 18
MAX COUNT       = 59
COUNTS OF 1     = 320412
OBSERVED COUNTS (60)
1	320412
2	3997
3	383
4	64
5	22
6	14
7	5
8	3
9	1
10	2
12	1
13	1
15	2
16	1
31	1
33	2
36	1
59	1

</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.45GB.
  PID: 419989;	Command: preseq;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_yield.txt']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_yield.txt`


> `preseq lc_extrap -v -o /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_yield.txt -B /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_dups_PE1.bam` (419991)

<pre>
BAM_INPUT
TOTAL READS     = 330356
DISTINCT READS  = 324913
DISTINCT COUNTS = 18
MAX COUNT       = 59
COUNTS OF 1     = 320412
MAX TERMS       = 10
OBSERVED COUNTS (60)
1	320412
2	3997
3	383
4	64
5	22
6	14
7	5
8	3
9	1
10	2
12	1
13	1
15	2
16	1
31	1
33	2
36	1
59	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
.._..._..._____._..._.____...._..___......_____.__..._.___...._.._.._...____.__.__..._..._.______..._.____._.......____._...__.._..___._._..__.._..._.._____.___..__.___._.____.__._._._..__..._______....
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.45GB.
  PID: 419991;	Command: preseq;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_counts.txt']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_counts.txt`


> `echo '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_yield.txt '$(samtools view -c -F 4 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_dups_PE1.bam)' '$(samtools view -c -F 4 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam) > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_counts.txt` (419994)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 419994;	Command: echo;	Return code: 0;	Memory used: 0.0GB

> `awk '{sum+=$2} END {printf "%.0f", sum}' /nv/t1/genomes/hg38/hg38.chrom.sizes`

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_plot.pdf', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_plot.png']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_plot.pdf`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_plot.png`


> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R preseq -i /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_yield.txt -r /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_counts.txt -o /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_preseq_plot` (420001)

<pre>
Processing tutorial
INFO: Found real counts for tutorial - Total: 330356 Unique: 329997

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.45GB.
  PID: 420001;	Command: Rscript;	Return code: 0;	Memory used: 0.0GB
> `Library complexity`	QC_hg38/tutorial_preseq_plot.pdf	Library complexity	QC_hg38/tutorial_preseq_plot.png	PEPPRO	_OBJ_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/objects.tsv

### Calculate NRF, PBC1, and PBC2 (07-31 20:54:46) elapsed: 11.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam.bai']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam.bai`


> `samtools index /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam` (420033)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420033;	Command: samtools;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_bamQC.tsv']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_bamQC.tsv`


> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py -i /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam -c 1 -o /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_bamQC.tsv` (420035)

<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam'
Temporary files will be stored in: '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tmp_tutorial_PE1_93hvn24v'
Processing with 1 cores...
[Name: chr1; Size: 248956422]
[Name: chr2; Size: 242193529]
[Name: chr3; Size: 198295559]
[Name: chr4; Size: 190214555]
[Name: chr5; Size: 181538259]
[Name: chr6; Size: 170805979]
[Name: chr7; Size: 159345973]
[Name: chr8; Size: 145138636]
[Name: chr9; Size: 138394717]
[Name: chr10; Size: 133797422]
[Name: chr11; Size: 135086622]
[Name: chr12; Size: 133275309]
[Name: chr13; Size: 114364328]
[Name: chr14; Size: 107043718]
[Name: chr15; Size: 101991189]
[Name: chr16; Size: 90338345]
[Name: chr17; Size: 83257441]
[Name: chr18; Size: 80373285]
[Name: chr19; Size: 58617616]
[Name: chr20; Size: 64444167]
[Name: chr21; Size: 46709983]
[Name: chr22; Size: 50818468]
[Name: chrX; Size: 156040895]
[Name: chrY; Size: 57227415]
[Name: chr1_KI270706v1_random; Size: 175055]
[Name: chr1_KI270713v1_random; Size: 40745]
[Name: chr14_GL000225v1_random; Size: 211173]
[Name: chr14_GL000194v1_random; Size: 191469]
[Name: chr14_KI270725v1_random; Size: 172810]
[Name: chr14_KI270726v1_random; Size: 43739]
[Name: chr17_GL000205v2_random; Size: 185591]
[Name: chr22_KI270733v1_random; Size: 179772]
[Name: chr22_KI270736v1_random; Size: 181920]
[Name: chrUn_KI270442v1; Size: 392061]
[Name: chrUn_KI270438v1; Size: 112505]
[Name: chrUn_GL000195v1; Size: 182896]
[Name: chrUn_GL000219v1; Size: 179198]
[Name: chrUn_KI270750v1; Size: 148850]
[Name: chrUn_GL000218v1; Size: 161147]
Discarding 156 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrEBV']
Keeping 39 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270713v1_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr17_GL000205v2_random', 'chr22_KI270733v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_KI270750v1', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 39 files into output file: '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_bamQC.tsv'
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.45GB.
  PID: 420035;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.0GB
Follow:

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_bamQC.tsv`


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_bamQC.tsv`


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_bamQC.tsv`


> `NRF`	1.0	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `PBC1`	164998.5	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `PBC2`	164998.5	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_unmap.bam']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_unmap.bam`


> `samtools view -b -@ 1 -f 12  /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_unmap.bam` (420082)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420082;	Command: samtools;	Return code: 0;	Memory used: 0.0GB
Follow:

> `samtools view -c -f 4 -@ 1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam`


> `Unmapped_reads`	44906	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

### Split BAM by strand (07-31 20:54:53) elapsed: 7.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam`


> `samtools view -bh -F 20 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam` (420086)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420086;	Command: samtools;	Return code: 0;	Memory used: 0.008GB

> `samtools view -bh -f 16 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam` (420089)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420089;	Command: samtools;	Return code: 0;	Memory used: 0.008GB

### Calculate TSS enrichment (07-31 20:54:55) elapsed: 3.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/plus_TSS.tsv', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/minus_TSS.tsv']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/plus_TSS.tsv`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/minus_TSS.tsv`


> `sed -n -e '/[[:space:]]+/w /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/minus_TSS.tsv' /nv/t1/genomes/hg38/refgene_anno/hg38_TSS.bed` (420092)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420092;	Command: sed;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_TssEnrichment.txt']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_TssEnrichment.txt`


> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/plus_TSS.tsv -p ends -c 1 -z -v -s 6 -o /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_TssEnrichment.txt` (420093)

<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.45GB.
  PID: 420093;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.066GB

> `TSS_Plus_Score`	56.2	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_TssEnrichment.txt']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_TssEnrichment.txt`


> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/minus_TSS.tsv -p ends -c 1 -z -v -s 6 -o /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_TssEnrichment.txt` (420142)

<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.45GB.
  PID: 420142;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.067GB

> `TSS_Minus_Score`	17.9	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_TSSenrichment.pdf']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_TSSenrichment.pdf`


> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_TssEnrichment.txt /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_TssEnrichment.txt` (420375)

<pre>

Generating TSS plot with /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_TssEnrichment.txt and /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.45GB.
  PID: 420375;	Command: Rscript;	Return code: 0;	Memory used: 0.0GB
> `TSS enrichment`	QC_hg38/tutorial_TSSenrichment.pdf	TSS enrichment	QC_hg38/tutorial_TSSenrichment.png	PEPPRO	_OBJ_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/objects.tsv
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt`


> `samtools view -H /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt` (420407,420408,420409,420410)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420408;	Command: grep;	Return code: 0;	Memory used: 0.0GB
  PID: 420407;	Command: samtools;	Return code: 0;	Memory used: 0.0GB
  PID: 420410;	Command: awk;	Return code: 0;	Memory used: 0.001GB
  PID: 420409;	Command: awk;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt`


> `cut -f 1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt` (420412)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420412;	Command: cut;	Return code: 0;	Memory used: 0.0GB

### Calculate Pause Index (PI) (07-31 20:55:11) elapsed: 16.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_tss.bed', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_gene_body.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_tss.bed`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_gene_body.bed`


> `grep -wf /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt /nv/t1/genomes/hg38/ensembl_gtf/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_tss.bed` (420414,420415)

<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.45GB.
  PID: 420414;	Command: grep;	Return code: 0;	Memory used: 0.0GB
  PID: 420415;	Command: bedtools;	Return code: 0;	Memory used: 0.098GB

> `grep -wf /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt /nv/t1/genomes/hg38/ensembl_gtf/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_gene_body.bed` (420419,420420)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420419;	Command: grep;	Return code: 0;	Memory used: 0.0GB
  PID: 420420;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_TSS_density.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_TSS_density.bed`


> `bedtools coverage -sorted -counts -s -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_tss.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_TSS_density.bed` (420423,420424,420425,420426)

<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.45GB.
  PID: 420423;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
  PID: 420425;	Command: sort;	Return code: 0;	Memory used: 0.0GB
  PID: 420424;	Command: awk;	Return code: 0;	Memory used: 0.0GB
  PID: 420426;	Command: sort;	Return code: 0;	Memory used: 0.007GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_gene_body_density.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_gene_body_density.bed`


> `bedtools coverage -sorted -counts -s -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_gene_body.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_gene_body_density.bed` (420429,420430,420431)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420430;	Command: awk;	Return code: 0;	Memory used: 0.0GB
  PID: 420429;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
  PID: 420431;	Command: sort;	Return code: 0;	Memory used: 0.001GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.txt']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.txt`


> `join -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_TSS_density.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_gene_body_density.bed | awk '{print ($6/($3-$2))/($9/($8-$7))}' > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.txt` (420434,420435)

join: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_gene_body_density.bed:33: is not sorted: chr16	15950077	16143074	ABCC1	.	+	73
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420434;	Command: join;	Return code: 1;	Memory used: 0.0GB
  PID: 420435;	Command: awk;	Return code: 0;	Memory used: 0.0GB
Subprocess returned nonzero result. Check above output for details
ERROR: Subprocess returned nonzero result, but pipeline is continuing because nofail=True

> `sort -n /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.txt | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`


> `Pause index`	108.18	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.pdf']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.pdf`


> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi -i /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.txt` (420440)

<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.45GB.
  PID: 420440;	Command: Rscript;	Return code: 0;	Memory used: 0.0GB
> `Pause index`	QC_hg38/tutorial_pause_index.pdf	Pause index	QC_hg38/tutorial_pause_index.png	PEPPRO	_OBJ_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/objects.tsv

### Calculate FRiP (07-31 20:55:22) elapsed: 11.0 _TIME_


> `samtools view -@ 4 -c -L /nv/t1/genomes/hg38/refgene_anno/hg38_pre-mRNA.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam`

326735.0 120762

> `Plus FRiP`	0.37	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `samtools view -@ 4 -c -L /nv/t1/genomes/hg38/refgene_anno/hg38_pre-mRNA.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam`

326735.0 116268

> `Minus FRiP`	0.36	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

### Plot fragment distribution (07-31 20:55:23) elapsed: 1.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLenDistribution.pdf']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLenDistribution.pdf`


> `perl /scratch/jps3dp/tools/databio/peppro/tools/fragment_length_dist.pl /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLen.txt` (420499)

<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.45GB.
  PID: 420499;	Command: perl;	Return code: 0;	Memory used: 0.003GB

> `sort -n  /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLen.txt | uniq -c  > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_fragCount.txt` (420502,420503)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420502;	Command: sort;	Return code: 0;	Memory used: 0.0GB
  PID: 420503;	Command: uniq;	Return code: 0;	Memory used: 0.0GB

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frag -l /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLen.txt -c /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_fragCount.txt -p /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLenDistribution.pdf -t /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLenDistribution.txt` (420506)

<pre>
Fragment distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.45GB.
  PID: 420506;	Command: Rscript;	Return code: 0;	Memory used: 0.0GB
> `Fragment distribution`	QC_hg38/tutorial_fragLenDistribution.pdf	Fragment distribution	QC_hg38/tutorial_fragLenDistribution.png	PEPPRO	_OBJ_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/objects.tsv
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed`


> `ln -sf /nv/t1/genomes/hg38/feat_annotation/hg38_annotations.bed.gz /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed.gz` (420536)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420536;	Command: ln;	Return code: 0;	Memory used: 0.0GB

> `gzip -f -d -c /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed.gz > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed` (420537)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420537;	Command: gzip;	Return code: 0;	Memory used: 0.0GB

### Calculate fraction of reads in features (FRiF) (07-31 20:55:30) elapsed: 7.0 _TIME_


> `cut -f 4 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed | sort -u`

original_path: ["/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3' UTR"]
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3' UTR`


> `awk -F'	' '{print>"/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/"$4}' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed` (420545)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420545;	Command: awk;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3_UTR']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3_UTR`


> `mv "/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3' UTR" "/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3_UTR"` (420547)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420547;	Command: mv;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3_UTR_sort.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3_UTR_sort.bed`


> `cut -f 1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | grep -wf - /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3_UTR_sort.bed` (420548,420549,420550,420551)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420548;	Command: cut;	Return code: 0;	Memory used: 0.0GB
  PID: 420549;	Command: grep;	Return code: 0;	Memory used: 0.0GB
  PID: 420551;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
  PID: 420550;	Command: cut;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_plus_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_plus_coverage.bed`


> `bedtools coverage -sorted -counts -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3_UTR_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_plus_coverage.bed` (420554)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420554;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_minus_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_minus_coverage.bed`


> `bedtools coverage -sorted -counts -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3_UTR_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_minus_coverage.bed` (420556)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420556;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
original_path: ["/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5' UTR"]
Target exists: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5' UTR`
Skipped command: `awk -F'	' '{print>"/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/"$4}' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed`
Command ID incremented by: `1`. Current ID: `117`

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5_UTR']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5_UTR`


> `mv "/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5' UTR" "/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5_UTR"` (420559)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420559;	Command: mv;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5_UTR_sort.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5_UTR_sort.bed`


> `cut -f 1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | grep -wf - /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5_UTR_sort.bed` (420560,420561,420562,420563)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420560;	Command: cut;	Return code: 0;	Memory used: 0.0GB
  PID: 420561;	Command: grep;	Return code: 0;	Memory used: 0.0GB
  PID: 420563;	Command: bedtools;	Return code: 0;	Memory used: 0.028GB
  PID: 420562;	Command: cut;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_plus_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_plus_coverage.bed`


> `bedtools coverage -sorted -counts -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5_UTR_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_plus_coverage.bed` (420565)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420565;	Command: bedtools;	Return code: 0;	Memory used: 0.001GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_minus_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_minus_coverage.bed`


> `bedtools coverage -sorted -counts -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5_UTR_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_minus_coverage.bed` (420568)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420568;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Enhancer']
Target exists: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Enhancer`
Skipped command: `awk -F'	' '{print>"/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/"$4}' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed`
Command ID incremented by: `1`. Current ID: `125`

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Enhancer_sort.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Enhancer_sort.bed`


> `cut -f 1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | grep -wf - /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Enhancer_sort.bed` (420570,420571,420572,420573)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420570;	Command: cut;	Return code: 0;	Memory used: 0.0GB
  PID: 420572;	Command: cut;	Return code: 0;	Memory used: 0.0GB
  PID: 420571;	Command: grep;	Return code: 0;	Memory used: 0.0GB
  PID: 420573;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_plus_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_plus_coverage.bed`


> `bedtools coverage -sorted -counts -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Enhancer_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_plus_coverage.bed` (420576)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420576;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_minus_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_minus_coverage.bed`


> `bedtools coverage -sorted -counts -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Enhancer_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_minus_coverage.bed` (420578)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420578;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Exon']
Target exists: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Exon`
Skipped command: `awk -F'	' '{print>"/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/"$4}' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed`
Command ID incremented by: `1`. Current ID: `132`

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Exon_sort.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Exon_sort.bed`


> `cut -f 1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | grep -wf - /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Exon_sort.bed` (420580,420581,420582,420583)

<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.45GB.
  PID: 420580;	Command: cut;	Return code: 0;	Memory used: 0.0GB
  PID: 420581;	Command: grep;	Return code: 0;	Memory used: 0.0GB
  PID: 420583;	Command: bedtools;	Return code: 0;	Memory used: 0.169GB
  PID: 420582;	Command: cut;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_plus_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_plus_coverage.bed`


> `bedtools coverage -sorted -counts -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Exon_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_plus_coverage.bed` (420588)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420588;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_minus_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_minus_coverage.bed`


> `bedtools coverage -sorted -counts -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Exon_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_minus_coverage.bed` (420591)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420591;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Intron']
Target exists: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Intron`
Skipped command: `awk -F'	' '{print>"/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/"$4}' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed`
Command ID incremented by: `1`. Current ID: `139`

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Intron_sort.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Intron_sort.bed`


> `cut -f 1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | grep -wf - /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Intron_sort.bed` (420593,420595,420596,420597)

<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.45GB.
  PID: 420593;	Command: cut;	Return code: 0;	Memory used: 0.0GB
  PID: 420597;	Command: bedtools;	Return code: 0;	Memory used: 0.081GB
  PID: 420595;	Command: grep;	Return code: 0;	Memory used: 0.0GB
  PID: 420596;	Command: cut;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_plus_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_plus_coverage.bed`


> `bedtools coverage -sorted -counts -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Intron_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_plus_coverage.bed` (420604)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420604;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_minus_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_minus_coverage.bed`


> `bedtools coverage -sorted -counts -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Intron_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_minus_coverage.bed` (420607)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420607;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter']
Target exists: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter`
Skipped command: `awk -F'	' '{print>"/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/"$4}' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed`
Command ID incremented by: `1`. Current ID: `146`

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_sort.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_sort.bed`


> `cut -f 1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | grep -wf - /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_sort.bed` (420609,420610,420611,420612)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420609;	Command: cut;	Return code: 0;	Memory used: 0.0GB
  PID: 420610;	Command: grep;	Return code: 0;	Memory used: 0.0GB
  PID: 420612;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB
  PID: 420611;	Command: cut;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_plus_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_plus_coverage.bed`


> `bedtools coverage -sorted -counts -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_plus_coverage.bed` (420615)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420615;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_minus_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_minus_coverage.bed`


> `bedtools coverage -sorted -counts -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_minus_coverage.bed` (420617)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420617;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter Flanking Region']
Target exists: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter Flanking Region`
Skipped command: `awk -F'	' '{print>"/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/"$4}' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed`
Command ID incremented by: `1`. Current ID: `153`

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region`


> `mv "/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter Flanking Region" "/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region"` (420619)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420619;	Command: mv;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region_sort.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region_sort.bed`


> `cut -f 1 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | grep -wf - /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region_sort.bed` (420620,420621,420622,420623)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420620;	Command: cut;	Return code: 0;	Memory used: 0.0GB
  PID: 420622;	Command: cut;	Return code: 0;	Memory used: 0.0GB
  PID: 420621;	Command: grep;	Return code: 0;	Memory used: 0.0GB
  PID: 420623;	Command: bedtools;	Return code: 0;	Memory used: 0.048GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_plus_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_plus_coverage.bed`


> `bedtools coverage -sorted -counts -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_plus_coverage.bed` (420626)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420626;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_minus_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_minus_coverage.bed`


> `bedtools coverage -sorted -counts -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_minus_coverage.bed` (420629)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420629;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB

### Plot FRiF (07-31 20:55:52) elapsed: 22.0 _TIME_


> `samtools view -@ 1 -q 10 -c -F4 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam`

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_frif.pdf']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_frif.pdf`


> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -n tutorial -r 166413 -o /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_frif.pdf --bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_plus_coverage.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_plus_coverage.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_plus_coverage.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_plus_coverage.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_plus_coverage.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_plus_coverage.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_plus_coverage.bed` (420632)

<pre>
Cumulative FRiF plot completed!

</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.45GB.
  PID: 420632;	Command: Rscript;	Return code: 0;	Memory used: 0.0GB
> `Plus FRiF`	QC_hg38/tutorial_plus_frif.pdf	Plus FRiF	QC_hg38/tutorial_plus_frif.png	PEPPRO	_OBJ_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/objects.tsv

> `samtools view -@ 1 -q 10 -c -F4 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam`

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_frif.pdf']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_frif.pdf`


> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -n tutorial -r 163584 -o /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_frif.pdf --bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_minus_coverage.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_minus_coverage.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_minus_coverage.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_minus_coverage.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_minus_coverage.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_minus_coverage.bed /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_minus_coverage.bed` (420681)

<pre>
Cumulative FRiF plot completed!

</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.45GB.
  PID: 420681;	Command: Rscript;	Return code: 0;	Memory used: 0.0GB
> `Minus FRiF`	QC_hg38/tutorial_minus_frif.pdf	Minus FRiF	QC_hg38/tutorial_minus_frif.png	PEPPRO	_OBJ_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/objects.tsv

### Calculate mRNA contamination (07-31 20:56:44) elapsed: 52.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_exons_sort.bed', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_introns_sort.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_exons_sort.bed`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_introns_sort.bed`


> `grep -wf /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt /nv/t1/genomes/hg38/refgene_anno/hg38_exons.bed | bedtools sort -i stdin -faidx /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_exons_sort.bed` (420730,420731)

<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.45GB.
  PID: 420731;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB
  PID: 420730;	Command: grep;	Return code: 0;	Memory used: 0.0GB

> `grep -wf /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt /nv/t1/genomes/hg38/refgene_anno/hg38_introns.bed | bedtools sort -i stdin -faidx /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_introns_sort.bed` (420737,420738,420739)

<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.45GB.
  PID: 420737;	Command: grep;	Return code: 0;	Memory used: 0.0GB
  PID: 420739;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB
  PID: 420738;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_coverage.bed', '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_coverage.bed']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_coverage.bed`,`/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_coverage.bed`


> `bedtools coverage -sorted -counts -s -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_exons_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_coverage.bed` (420747)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420747;	Command: bedtools;	Return code: 0;	Memory used: 0.0GB

> `bedtools coverage -sorted -counts -s -a /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_introns_sort.bed -b /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_PE1.bam -g /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_coverage.bed` (420750)

<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.45GB.
  PID: 420750;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_rpkm.tsv']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_rpkm.tsv`


> `awk -v OFS='	' '{readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4} END { for (a in readCount) { print gene[a], (readCount[a]/0.326735)/geneSizeKB[a]}}' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_coverage.bed | awk '$2>0' | sort -k1 > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_rpkm.tsv` (420753,420754,420755)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420753;	Command: awk;	Return code: 0;	Memory used: 0.0GB
  PID: 420755;	Command: sort;	Return code: 0;	Memory used: 0.0GB
  PID: 420754;	Command: awk;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_rpkm.tsv']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_rpkm.tsv`


> `awk -v OFS='	' '{readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4} END { for (a in readCount) { print gene[a], (readCount[a]/0.326735)/geneSizeKB[a]}}' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_coverage.bed | awk '$2>0' | sort -k1 > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_rpkm.tsv` (420757,420758,420759)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420757;	Command: awk;	Return code: 0;	Memory used: 0.0GB
  PID: 420759;	Command: sort;	Return code: 0;	Memory used: 0.0GB
  PID: 420758;	Command: awk;	Return code: 0;	Memory used: 0.0GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_intron_exon.tsv']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_intron_exon.tsv`


> `join -a1 -a2 -j1 -e0 -o 0 1.2 2.2 /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_rpkm.tsv /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_rpkm.tsv > /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_intron_exon.tsv` (420761)

<pre>
join: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_rpkm.tsv:12: is not sorted: AARS	2.62812
join: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_rpkm.tsv:8: is not sorted: AARS	2.82081
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420761;	Command: join;	Return code: 1;	Memory used: 0.0GB
Subprocess returned nonzero result. Check above output for details
ERROR: Subprocess returned nonzero result, but pipeline is continuing because nofail=True

> `awk '$2>0' /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_intron_exon.tsv | awk '{print $3/$2}' | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`


> `mRNA contamination`	0.63	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_mRNA_contamination.pdf']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_mRNA_contamination.pdf`


> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_intron_exon.tsv --raw` (420768)

<pre>
mRNA contamination plot completed!

Warning messages:
1: Removed 3008 rows containing non-finite values (stat_boxplot). 
2: Removed 3008 rows containing non-finite values (stat_boxplot). 
3: Removed 3008 rows containing non-finite values (stat_summary). 
4: Removed 3008 rows containing non-finite values (stat_boxplot). 
5: Removed 3008 rows containing non-finite values (stat_boxplot). 
6: Removed 3008 rows containing non-finite values (stat_summary). 
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.45GB.
  PID: 420768;	Command: Rscript;	Return code: 0;	Memory used: 0.0GB
> `mRNA contamination`	QC_hg38/tutorial_mRNA_contamination.pdf	mRNA contamination	QC_hg38/tutorial_mRNA_contamination.png	PEPPRO	_OBJ_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/objects.tsv

### Produce bigWig files (07-31 20:57:04) elapsed: 20.0 _TIME_

original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/signal_hg38/tutorial_plus_body_0-mer.bw']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/signal_hg38/tutorial_plus_body_0-mer.bw`


> `samtools index /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam` (420816)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 420816;	Command: samtools;	Return code: 0;	Memory used: 0.0GB

> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam -c /nv/t1/genomes/hg38/hg38.chrom.sizes -o /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/signal_hg38/tutorial_plus_body_0-mer.bw -p 1 --variable-step --tail-edge` (420817)

<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_plus.bam'
Temporary files will be stored in: 'tmp_tutorial_plus_v2s7uhd7'
Processing with 1 cores...
[Name: chr1; Size: 248956422]
[Name: chr2; Size: 242193529]
[Name: chr3; Size: 198295559]
[Name: chr4; Size: 190214555]
[Name: chr5; Size: 181538259]
[Name: chr6; Size: 170805979]
[Name: chr7; Size: 159345973]
[Name: chr8; Size: 145138636]
[Name: chr9; Size: 138394717]
[Name: chr10; Size: 133797422]
[Name: chr11; Size: 135086622]
[Name: chr12; Size: 133275309]
[Name: chr13; Size: 114364328]
[Name: chr14; Size: 107043718]
[Name: chr15; Size: 101991189]
[Name: chr16; Size: 90338345]
[Name: chr17; Size: 83257441]
[Name: chr18; Size: 80373285]
[Name: chr19; Size: 58617616]
[Name: chr20; Size: 64444167]
[Name: chr21; Size: 46709983]
[Name: chr22; Size: 50818468]
[Name: chrX; Size: 156040895]
[Name: chrY; Size: 57227415]
[Name: chr1_KI270706v1_random; Size: 175055]
[Name: chr1_KI270713v1_random; Size: 40745]
[Name: chr14_KI270726v1_random; Size: 43739]
[Name: chr17_GL000205v2_random; Size: 185591]
[Name: chr22_KI270733v1_random; Size: 179772]
[Name: chr22_KI270736v1_random; Size: 181920]
[Name: chrUn_KI270442v1; Size: 392061]
[Name: chrUn_GL000195v1; Size: 182896]
[Name: chrUn_GL000219v1; Size: 179198]
[Name: chrUn_KI270750v1; Size: 148850]
Discarding 161 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Keeping 34 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270713v1_random', 'chr14_KI270726v1_random', 'chr17_GL000205v2_random', 'chr22_KI270733v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_KI270750v1']
Reduce step (merge files)...
Merging 34 files into output file: '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/signal_hg38/tutorial_plus_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.45GB.
  PID: 420817;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.006GB
original_path: ['/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/signal_hg38/tutorial_minus_body_0-mer.bw']
Target to produce: `/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/signal_hg38/tutorial_minus_body_0-mer.bw`


> `samtools index /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam` (421014)

<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.45GB.
  PID: 421014;	Command: samtools;	Return code: 0;	Memory used: 0.0GB

> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam -c /nv/t1/genomes/hg38/hg38.chrom.sizes -o /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/signal_hg38/tutorial_minus_body_0-mer.bw -p 1 --variable-step --tail-edge` (421015)

<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_minus.bam'
Temporary files will be stored in: 'tmp_tutorial_minus_lycz3ftu'
Processing with 1 cores...
[Name: chr1; Size: 248956422]
[Name: chr2; Size: 242193529]
[Name: chr3; Size: 198295559]
[Name: chr4; Size: 190214555]
[Name: chr5; Size: 181538259]
[Name: chr6; Size: 170805979]
[Name: chr7; Size: 159345973]
[Name: chr8; Size: 145138636]
[Name: chr9; Size: 138394717]
[Name: chr10; Size: 133797422]
[Name: chr11; Size: 135086622]
[Name: chr12; Size: 133275309]
[Name: chr13; Size: 114364328]
[Name: chr14; Size: 107043718]
[Name: chr15; Size: 101991189]
[Name: chr16; Size: 90338345]
[Name: chr17; Size: 83257441]
[Name: chr18; Size: 80373285]
[Name: chr19; Size: 58617616]
[Name: chr20; Size: 64444167]
[Name: chr21; Size: 46709983]
[Name: chr22; Size: 50818468]
[Name: chrX; Size: 156040895]
[Name: chrY; Size: 57227415]
[Name: chr1_KI270713v1_random; Size: 40745]
[Name: chr14_GL000225v1_random; Size: 211173]
[Name: chr14_GL000194v1_random; Size: 191469]
[Name: chr14_KI270725v1_random; Size: 172810]
[Name: chr14_KI270726v1_random; Size: 43739]
[Name: chr17_GL000205v2_random; Size: 185591]
[Name: chr22_KI270733v1_random; Size: 179772]
[Name: chrUn_KI270438v1; Size: 112505]
[Name: chrUn_GL000195v1; Size: 182896]
[Name: chrUn_GL000219v1; Size: 179198]
[Name: chrUn_KI270750v1; Size: 148850]
[Name: chrUn_GL000218v1; Size: 161147]
Discarding 159 chunk(s) of reads: ['chrM', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrEBV']
Keeping 36 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270713v1_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr17_GL000205v2_random', 'chr22_KI270733v1_random', 'chrUn_KI270438v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_KI270750v1', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 36 files into output file: '/sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/signal_hg38/tutorial_minus_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.45GB.
  PID: 421015;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.007GB

Changed status from running to completed.
Starting cleanup: 79 files; 8 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_noadap.fastq
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_noadap.fastq`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_dedup.fastq
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_dedup.fastq`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq.single.fq
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq.single.fq`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq.single.fq
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq.single.fq`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/repaired.flag

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq.single.fq
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_trimmed.fastq.single.fq`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq.single.fq
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq.single.fq`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/dups_repaired.flag

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tmpcmta8okh
`rmdir /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tmpcmta8okh`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/human_rDNA_bt2

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tmp4z58ml5o
`rmdir /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tmp4z58ml5o`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/human_rDNA_dups_bt2

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tmp8zeo8d45
`rmdir /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tmp8zeo8d45`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/rCRSd_bt2

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tmpd16_9pvj
`rmdir /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tmpd16_9pvj`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/rCRSd_dups_bt2

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tmplh2_a83a
`rmdir /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tmplh2_a83a`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_fail_qc_dups.bam
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_fail_qc_dups.bam`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_dups_R1.fq
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_dups_R1.fq`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_dups_R2.fq
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_human_rDNA_unmap_dups_R2.fq`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_dups_R1.fq
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_dups_R1.fq`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_dups_R2.fq
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/tutorial_rCRSd_unmap_dups_R2.fq`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam.bai
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam.bai`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam.bai
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort.bam.bai`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp_dups.bam.bai
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp_dups.bam.bai`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort_dups.bam.bai
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_sort_dups.bam.bai`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_dups_PE1.bam
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_dups_PE1.bam`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp_dups.bam
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp_dups.bam`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/aligned_hg38/tutorial_temp.bam`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/plus_TSS.tsv
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/plus_TSS.tsv`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_TssEnrichment.txt
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_plus_TssEnrichment.txt`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/minus_TSS.tsv
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/minus_TSS.tsv`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_TssEnrichment.txt
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_minus_TssEnrichment.txt`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_order.txt`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/chr_keep.txt`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_tss.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_tss.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_gene_body.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_ensembl_gene_body.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_TSS_density.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_TSS_density.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_gene_body_density.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_gene_body_density.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.txt
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_pause_index.txt`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLen.txt
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_fragLen.txt`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_fragCount.txt
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_fragCount.txt`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/raw/hg38_annotations.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3_UTR
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3_UTR`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3_UTR_sort.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/3_UTR_sort.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_plus_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_plus_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_minus_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_3_UTR_minus_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5_UTR
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5_UTR`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5_UTR_sort.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/5_UTR_sort.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_plus_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_plus_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_minus_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_5_UTR_minus_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Enhancer
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Enhancer`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Enhancer_sort.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Enhancer_sort.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_plus_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_plus_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_minus_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Enhancer_minus_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Exon
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Exon`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Exon_sort.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Exon_sort.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_plus_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_plus_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_minus_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Exon_minus_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Intron
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Intron`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Intron_sort.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Intron_sort.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_plus_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_plus_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_minus_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Intron_minus_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_sort.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_sort.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_plus_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_plus_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_minus_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_minus_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region_sort.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/Promoter_Flanking_Region_sort.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_plus_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_plus_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_minus_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_Promoter_Flanking_Region_minus_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_exons_sort.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_exons_sort.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_introns_sort.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/hg38_introns_sort.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_coverage.bed
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_coverage.bed`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_rpkm.tsv
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_exons_rpkm.tsv`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_rpkm.tsv
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_introns_rpkm.tsv`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_intron_exon.tsv
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/QC_hg38/tutorial_intron_exon.tsv`

Cleaning up conditional list. . .

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial*.fastq
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1_processed.fastq`
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R1.fastq`
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2.fastq`
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed_dups.fastq`
`rm /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/tutorial_R2_trimmed.fastq`

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/*.fq

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/*.fastq

Removing glob: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/fastq/*.log

Removing glob: mkfifo /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/human_rDNA_bt2

Removing glob: mkfifo /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/human_rDNA_dups_bt2

Removing glob: mkfifo /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/rCRSd_bt2

Removing glob: mkfifo /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/prealignments/rCRSd_dups_bt2

> `Time`	0:14:01	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

> `Success`	07-31-20:57:16	PEPPRO	_RES_
original_path: /sfs/lustre/allocations/shefflab/processed/peppro_tutorial/pe/07-31-19/results_pipeline/tutorial/stats.tsv

##### [Epilogue:]
*   Total elapsed time:  0:14:01
*     Peak memory used:  3.45 GB
* Pipeline completed at:  (07-31 20:57:16) elapsed: 13.0 _TIME_
