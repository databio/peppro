# PEPPRO Test Suite

This directory contains the PEPPRO test suite, organized into two tiers:

- **Unit tests** — fast, no genome data or external bioinformatics tools required; run on every push/PR via GitHub Actions
- **Integration tests** — full pipeline runs; require a self-hosted runner with genome indices and all tools installed

---

## Directory Structure

```
tests/
├── data/                       # Small test FASTQ files (~3 MB total)
│   ├── test_R1.fastq.gz        # SE reads (12,500 reads)
│   ├── test_R2.fastq.gz        # PE reverse reads (rev-comp of R1)
│   └── test_R1_umi.fastq.gz    # R1 with 8-nt UMI prefix for UMI tests
├── pep_configs/                # PEP project configs for each scenario
│   ├── se_basic.yaml / .csv
│   ├── pe_basic.yaml / .csv
│   └── ...
├── looper_configs/             # Looper run configs for each scenario
│   ├── .looper_se_basic.yaml
│   └── ...
├── scripts/
│   └── generate_test_data.sh   # Regenerate test FASTQ data from source
├── test_unit.py                # Unit tests (no tools/genome needed)
├── test_integration.py         # Integration tests (full pipeline runs)
└── README.md                   # This file
```

---

## Unit Tests

Unit tests cover:

- **Constants**: `RUNON_SOURCE`, `ADAPTER_REMOVERS`, `TRIMMERS`, `DEDUPLICATORS` values and defaults
- **PEP loading**: Each test config loads correctly with expected sample attributes
- **Schema validation**: eido validation passes for valid configs; regression tests ensure invalid inputs (e.g., integer `umi_len` in YAML `imply`, invalid `protocol`/`adapter`/`trimmer`/`dedup` enum values) fail correctly
- **Argument parsing**: All CLI flags parse correctly, defaults are correct, invalid choices raise `SystemExit`
- **Recovery paths**: Expected output file naming conventions are documented and verified

### Running unit tests

```bash
# Via pytest directly
pytest tests/test_unit.py -v

# Via Makefile
make test-unit
```

No environment variables or external tools are needed.

---

## Integration Tests

Integration tests run the full PEPPRO pipeline for each scenario and verify:

1. Pipeline exits with status `0`
2. Key output files exist (BAM, bigWig, stats.yaml)
3. `stats.yaml` contains the expected result keys
4. The `TestRecovery` class additionally tests checkpoint skipping and the `unmap_R1.fq` recovery regression

### Prerequisites

The integration tests require a machine with all PEPPRO dependencies installed and genome assets configured via refgenie:

| Tool | Version tested |
|------|---------------|
| bowtie2 | ≥2.4 |
| samtools | ≥1.13 |
| bedtools | ≥2.30 |
| cutadapt | ≥4.0 |
| fastp | ≥0.23 |
| seqtk | ≥1.3 |
| fastx_toolkit | any |
| seqkit | ≥2.0 |
| fqdedup | any |
| fastq_pair | any |
| wigToBigWig | UCSC |
| bedGraphToBigWig | UCSC |

**Genome assets** (via refgenie, pointed to by `$REFGENIE`):

- `hg38/bowtie2_index`
- `human_rDNA/bowtie2_index`
- `hg38/fasta` (for chromosome sizes)
- `hg38/blacklist` (optional, for coverage tests)

### Running integration tests

**Important notes:**

- The PyPI package for pypiper is **`piper`** (not `pypiper`, which is an unrelated package).
- Bioinformatics tools (samtools, bowtie2, etc.) are provided via bulker. The wrapper script handles this automatically, or you can use `bulker activate` / `bulker exec` directly.
- Tests run with `-p local` (divvy local compute package) so the pipeline executes inline rather than being submitted to a job scheduler.

```bash
# Recommended: use the wrapper script (runs pytest via bulker exec)
bash tests/scripts/test-integration.sh

# Or manually: activate bulker, then run pytest
bulker activate databio/peppro:1.1.0
RUN_INTEGRATION_TESTS=true pytest tests/test_integration.py -v
bulker deactivate

# Run a specific scenario
bash tests/scripts/test-integration.sh -k se_basic

# Via Makefile targets
make test-se          # All SE scenarios
make test-pe          # All PE scenarios
make test-recovery    # Recovery regression tests
make test-integration # All integration tests
make test-all         # Unit + integration

# Run a single named scenario
make test-scenario SCENARIO=se_fastp

# Keep output directories for debugging (default: cleaned up after each class)
KEEP_TEST_OUTPUTS=true RUN_INTEGRATION_TESTS=true pytest tests/test_integration.py -v -k se_basic
```

---

## Test Scenarios

| Scenario | Read type | Protocol | Adapter | Trimmer | Dedup | Notes |
|----------|-----------|----------|---------|---------|-------|-------|
| `se_basic` | SE | PRO-seq | cutadapt | seqtk | — | Baseline SE run |
| `pe_basic` | PE | PRO-seq | cutadapt | seqtk | — | Baseline PE run |
| `se_groseq` | SE | GRO-seq | cutadapt | seqtk | — | GRO-seq protocol |
| `se_umi` | SE | PRO-seq | cutadapt | seqtk | seqkit | 8-nt UMI deduplication |
| `pe_umi` | PE | PRO-seq | cutadapt | seqtk | seqkit | PE with UMI dedup |
| `se_fastp` | SE | PRO-seq | fastp | seqtk | — | fastp adapter trimming |
| `se_fastx` | SE | PRO-seq | cutadapt | fastx | — | fastx_trimmer |
| `se_fqdedup` | SE | PRO-seq | cutadapt | seqtk | fqdedup | fqdedup UMI dedup |
| `se_scale` | SE | PRO-seq | cutadapt | seqtk | — | `--scale` flag |
| `se_no_complexity` | SE | PRO-seq | cutadapt | seqtk | — | `--no-complexity` flag |
| `se_nofifo` | SE | PRO-seq | cutadapt | seqtk | — | `--no-fifo` flag |
| `se_coverage` | SE | PRO-seq | cutadapt | seqtk | — | `--coverage` flag |

---

## Test Data

The files in `tests/data/` are derived from `examples/data/test_r1.fq.gz` (the existing pipeline example read file). They are small enough to commit to the repository (~1 MB each).

To regenerate the test data files (requires `seqtk`):

```bash
make test-data
# or
bash tests/scripts/generate_test_data.sh
```

---

## GitHub Actions

Unit tests run automatically on every push and pull request targeting `master` or `dev`, across Python 3.9, 3.11, and 3.12.

Integration tests are triggered manually via **workflow_dispatch** on a self-hosted runner:

1. Go to **Actions** → **Tests** → **Run workflow**
2. Set "Run integration tests" to `true`
3. Click **Run workflow**

See `.github/workflows/tests.yml` for the full configuration.
