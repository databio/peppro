test:
	python pipelines/peppro.py  -P 3 -M 100 -O test_out -R -S liver -G hg19  -Q paired  -C peppro.yaml  --genome-size hs --prealignments rCRSd human_repeats -I examples/test_data/liver-CD31_test_R1.fastq.gz -I2 examples/test_data/liver-CD31_test_R2.fastq.gz  

docker:
	docker build -t databio/peppro -f containers/peppro.Dockerfile .

singularity:
	singularity build ${SIMAGES}peppro docker://databio/peppro