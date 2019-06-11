test:
	python pipelines/peppro.py  -P 3 -M 100 -O peppro_test -R -S test -G hg38  -Q single  -C peppro.yaml  --genome-size hs --prealignments rCRSd human_repeats -I examples/data/test_R1.fq.gz

docker:
	docker build -t databio/peppro -f containers/peppro.Dockerfile .

singularity:
	singularity build $${SIMAGES}peppro docker://databio/peppro