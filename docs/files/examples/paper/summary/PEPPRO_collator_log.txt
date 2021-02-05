### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro_collator.py --config /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/peppro_paper_rerun.yaml -O /project/shefflab/processed//peppro/paper/6.11.2020 -P 1 -M 6000 -n PEPPRO -r /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         Compute host:  udc-ba34-36
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/summary/
*  Pipeline started at:   (06-15 13:09:28) elapsed: 14.0 _TIME_

### Version log:

*       Python version:  3.6.6
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.0.3
*        Pipeline hash:  887a9a85408962dc3ca070dc954e59a5d4d73a10
*      Pipeline branch:  * master
*        Pipeline date:  2020-06-09 11:36:49 -0400
*        Pipeline diff:  40 files changed, 16413 insertions(+), 3702 deletions(-)

### Arguments passed to pipeline:

*        `config_file`:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/peppro_paper_rerun.yaml`
*              `cores`:  `1`
*              `dirty`:  `False`
*       `force_follow`:  `False`
*             `logdev`:  `False`
*                `mem`:  `6000`
*               `name`:  `PEPPRO`
*          `new_start`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020`
*            `recover`:  `False`
*            `results`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*             `silent`:  `False`
*           `testmode`:  `False`
*          `verbosity`:  `None`

----------------------------------------

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/summary/PEPPRO_libComplexity.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/summary/PEPPRO_countData.csv`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO_summarizer.R /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/peppro_paper_rerun.yaml /project/shefflab/processed//peppro/paper/6.11.2020 /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline` (124768)
<pre>
Loading config file: /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/peppro_paper_rerun.yaml
Creating assets summary...
Summary (n=47): /project/shefflab/processed//peppro/paper/6.11.2020/PEPPRO_assets_summary.tsv
17 of 47 library complexity files available.
Processing Jurkat_ChRO-seq_1
Processing Jurkat_ChRO-seq_2
Processing HEK_PRO-seq
Processing HEK_ARF_PRO-seq
Processing H9_treated_PRO-seq_1
Processing H9_treated_PRO-seq_2
Processing H9_treated_PRO-seq_3
Processing H9_PRO-seq_10
Processing H9_PRO-seq_20
Processing H9_PRO-seq_30
Processing H9_PRO-seq_40
Processing H9_PRO-seq_50
Processing H9_PRO-seq_60
Processing H9_PRO-seq_70
Processing H9_PRO-seq_80
Processing H9_PRO-seq_90
Processing H9_PRO-seq_100
INFO: Found real counts for Jurkat_ChRO-seq_1 - Total (M): 21.334642 Unique (M): 13.139331
INFO: Found real counts for Jurkat_ChRO-seq_2 - Total (M): 33.989659 Unique (M): 30.4684
INFO: Found real counts for HEK_PRO-seq - Total (M): 15.183469 Unique (M): 14.923343
INFO: Found real counts for HEK_ARF_PRO-seq - Total (M): 10.869678 Unique (M): 10.681803
INFO: Found real counts for H9_treated_PRO-seq_1 - Total (M): 19.73867 Unique (M): 18.664481
INFO: Found real counts for H9_treated_PRO-seq_2 - Total (M): 20.247016 Unique (M): 19.595024
INFO: Found real counts for H9_treated_PRO-seq_3 - Total (M): 22.825749 Unique (M): 20.168929
INFO: Found real counts for H9_PRO-seq_10 - Total (M): 1.945718 Unique (M): 1.966383
INFO: Found real counts for H9_PRO-seq_20 - Total (M): 3.896875 Unique (M): 3.908887
INFO: Found real counts for H9_PRO-seq_30 - Total (M): 5.844297 Unique (M): 5.818689
INFO: Found real counts for H9_PRO-seq_40 - Total (M): 7.79292 Unique (M): 7.703559
INFO: Found real counts for H9_PRO-seq_50 - Total (M): 9.741802 Unique (M): 9.562095
INFO: Found real counts for H9_PRO-seq_60 - Total (M): 11.690255 Unique (M): 11.393512
INFO: Found real counts for H9_PRO-seq_70 - Total (M): 13.641283 Unique (M): 13.202733
INFO: Found real counts for H9_PRO-seq_80 - Total (M): 15.588108 Unique (M): 14.983475
INFO: Found real counts for H9_PRO-seq_90 - Total (M): 17.534903 Unique (M): 16.739307
INFO: Found real counts for H9_PRO-seq_100 - Total (M): 19.4822 Unique (M): 18.471612

43 of 47 gene counts files available
Counts table: /project/shefflab/processed//peppro/paper/6.11.2020/summary/PEPPRO_countData.csv

</pre>
Command completed. Elapsed time: 0:00:45. Running peak memory: 0.526GB.  
  PID: 124768;	Command: Rscript;	Return code: 0;	Memory used: 0.526GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:00:59
*  Total elapsed time (all runs):  0:00:45
*         Peak memory (this run):  0.5263 GB
*        Pipeline completed time: 2020-06-15 13:10:13
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro_collator.py --config /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/peppro_paper_rerun.yaml -O /project/shefflab/processed//peppro/paper/6.11.2020 -P 1 -M 6000 -n PEPPRO -r /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         Compute host:  udc-ba34-36
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/summary/
*  Pipeline started at:   (06-16 09:04:50) elapsed: 31.0 _TIME_

### Version log:

*       Python version:  3.6.6
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.0.3
*        Pipeline hash:  887a9a85408962dc3ca070dc954e59a5d4d73a10
*      Pipeline branch:  * master
*        Pipeline date:  2020-06-09 11:36:49 -0400
*        Pipeline diff:  40 files changed, 16413 insertions(+), 3702 deletions(-)

### Arguments passed to pipeline:

*        `config_file`:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/peppro_paper_rerun.yaml`
*              `cores`:  `1`
*              `dirty`:  `False`
*       `force_follow`:  `False`
*             `logdev`:  `False`
*                `mem`:  `6000`
*               `name`:  `PEPPRO`
*          `new_start`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020`
*            `recover`:  `False`
*            `results`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*             `silent`:  `False`
*           `testmode`:  `False`
*          `verbosity`:  `None`

----------------------------------------

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/summary/PEPPRO_libComplexity.pdf`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/summary/PEPPRO_countData.csv`  

### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:00:31
*  Total elapsed time (all runs):  0:00:45
*         Peak memory (this run):  0 GB
*        Pipeline completed time: 2020-06-16 09:04:50
Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/summary/PEPPRO collator_completed.flag'
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro_collator.py --config /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/peppro_paper_rerun.yaml -O /project/shefflab/processed//peppro/paper/6.11.2020 -P 1 -M 6000 -n PEPPRO -r /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline --new-start`
*         Compute host:  udc-ba34-36
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/summary/
*  Pipeline started at:   (06-16 09:05:10) elapsed: 0.0 _TIME_

### Version log:

*       Python version:  3.6.6
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.0.3
*        Pipeline hash:  887a9a85408962dc3ca070dc954e59a5d4d73a10
*      Pipeline branch:  * master
*        Pipeline date:  2020-06-09 11:36:49 -0400
*        Pipeline diff:  40 files changed, 16413 insertions(+), 3702 deletions(-)

### Arguments passed to pipeline:

*        `config_file`:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/peppro_paper_rerun.yaml`
*              `cores`:  `1`
*              `dirty`:  `False`
*       `force_follow`:  `False`
*             `logdev`:  `False`
*                `mem`:  `6000`
*               `name`:  `PEPPRO`
*          `new_start`:  `True`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020`
*            `recover`:  `False`
*            `results`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*             `silent`:  `False`
*           `testmode`:  `False`
*          `verbosity`:  `None`

----------------------------------------

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/summary/PEPPRO_libComplexity.pdf`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/summary/PEPPRO_countData.csv`  
New start mode; run anyway.  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/summary/PEPPRO_libComplexity.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/summary/PEPPRO_countData.csv`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO_summarizer.R /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/peppro_paper_rerun.yaml /project/shefflab/processed//peppro/paper/6.11.2020 /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline --new-start` (132188)
<pre>
Loading config file: /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/peppro_paper_rerun.yaml
Creating assets summary...
Summary (n=47): /project/shefflab/processed//peppro/paper/6.11.2020/PEPPRO_assets_summary.tsv
20 of 47 library complexity files available.
Processing Jurkat_ChRO-seq_1
Processing Jurkat_ChRO-seq_2
Processing HEK_PRO-seq
Processing HEK_ARF_PRO-seq
Processing H9_PRO-seq_1
Processing H9_PRO-seq_2
Processing H9_PRO-seq_3
Processing H9_treated_PRO-seq_1
Processing H9_treated_PRO-seq_2
Processing H9_treated_PRO-seq_3
Processing H9_PRO-seq_10
Processing H9_PRO-seq_20
Processing H9_PRO-seq_30
Processing H9_PRO-seq_40
Processing H9_PRO-seq_50
Processing H9_PRO-seq_60
Processing H9_PRO-seq_70
Processing H9_PRO-seq_80
Processing H9_PRO-seq_90
Processing H9_PRO-seq_100
INFO: Found real counts for Jurkat_ChRO-seq_1 - Total (M): 21.334642 Unique (M): 13.139331
INFO: Found real counts for Jurkat_ChRO-seq_2 - Total (M): 33.989659 Unique (M): 30.4684
INFO: Found real counts for HEK_PRO-seq - Total (M): 15.183469 Unique (M): 14.923343
INFO: Found real counts for HEK_ARF_PRO-seq - Total (M): 10.869678 Unique (M): 10.681803
INFO: Found real counts for H9_PRO-seq_1 - Total (M): 21.966708 Unique (M): 20.33023
INFO: Found real counts for H9_PRO-seq_2 - Total (M): 19.4822 Unique (M): 18.471612
INFO: Found real counts for H9_PRO-seq_3 - Total (M): 19.481041 Unique (M): 19.005606
INFO: Found real counts for H9_treated_PRO-seq_1 - Total (M): 19.73867 Unique (M): 18.664481
INFO: Found real counts for H9_treated_PRO-seq_2 - Total (M): 20.247016 Unique (M): 19.595024
INFO: Found real counts for H9_treated_PRO-seq_3 - Total (M): 22.825749 Unique (M): 20.168929
INFO: Found real counts for H9_PRO-seq_10 - Total (M): 1.945718 Unique (M): 1.966383
INFO: Found real counts for H9_PRO-seq_20 - Total (M): 3.896875 Unique (M): 3.908887
INFO: Found real counts for H9_PRO-seq_30 - Total (M): 5.844297 Unique (M): 5.818689
INFO: Found real counts for H9_PRO-seq_40 - Total (M): 7.79292 Unique (M): 7.703559
INFO: Found real counts for H9_PRO-seq_50 - Total (M): 9.741802 Unique (M): 9.562095
INFO: Found real counts for H9_PRO-seq_60 - Total (M): 11.690255 Unique (M): 11.393512
INFO: Found real counts for H9_PRO-seq_70 - Total (M): 13.641283 Unique (M): 13.202733
INFO: Found real counts for H9_PRO-seq_80 - Total (M): 15.588108 Unique (M): 14.983475
INFO: Found real counts for H9_PRO-seq_90 - Total (M): 17.534903 Unique (M): 16.739307
INFO: Found real counts for H9_PRO-seq_100 - Total (M): 19.4822 Unique (M): 18.471612

47 of 47 gene counts files available
Counts table: /project/shefflab/processed//peppro/paper/6.11.2020/summary/PEPPRO_countData.csv

</pre>
Command completed. Elapsed time: 0:00:48. Running peak memory: 0.567GB.  
  PID: 132188;	Command: Rscript;	Return code: 0;	Memory used: 0.567GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:00:48
*  Total elapsed time (all runs):  0:00:48
*         Peak memory (this run):  0.5671 GB
*        Pipeline completed time: 2020-06-16 09:05:58
Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/summary/PEPPRO collator_completed.flag'
