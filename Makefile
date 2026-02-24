test:
	python pipelines/peppro.py -P 3 -M 100 -O peppro_test -R -S test -G hg38 \
		-Q single -C peppro.yaml \
		--protocol PRO \
		--prealignment-names rCRSd human_repeats \
		--genome-index $$(refgenie seek hg38/bowtie2_index --seek-key dir) \
		--chrom-sizes $$(refgenie seek hg38/fasta --seek-key chrom_sizes) \
		--pipestat-schema peppro_output_schema.yaml \
		-I examples/data/test_r1.fq.gz

# -----------------------------------------------------------------------
# Test suite targets
# -----------------------------------------------------------------------

# Run only unit tests (no genome data or external tools required)
test-unit:
	pytest tests/test_unit.py -v --tb=short

# Run a single integration scenario (e.g. make test-se SCENARIO=se_basic)
SCENARIO ?= se_basic
test-scenario:
	RUN_INTEGRATION_TESTS=true pytest tests/test_integration.py -v --tb=short -k "$(SCENARIO)"

# Run all SE integration scenarios
test-se:
	RUN_INTEGRATION_TESTS=true pytest tests/test_integration.py -v --tb=short \
		-k "se_basic or se_groseq or se_umi or se_fastp or se_fastx or se_fqdedup or se_scale or se_no_complexity or se_nofifo or se_coverage"

# Run all PE integration scenarios
test-pe:
	RUN_INTEGRATION_TESTS=true pytest tests/test_integration.py -v --tb=short \
		-k "pe_basic or pe_umi"

# Run recovery regression tests
test-recovery:
	RUN_INTEGRATION_TESTS=true pytest tests/test_integration.py -v --tb=short \
		-k "recovery"

# Run all integration tests (SE + PE + recovery)
test-integration:
	RUN_INTEGRATION_TESTS=true pytest tests/test_integration.py -v --tb=short

# Run both unit and integration tests
test-all:
	RUN_INTEGRATION_TESTS=true pytest tests/ -v --tb=short

# Regenerate test FASTQ data files from the source R1 read file
test-data:
	bash tests/scripts/generate_test_data.sh

docker:
	docker build -t databio/peppro -f containers/peppro.Dockerfile .

singularity:
	singularity build $${SIMAGES}peppro docker://databio/peppro