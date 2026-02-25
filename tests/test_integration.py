"""
Integration tests for PEPPRO — runs the full pipeline for each test scenario.

Prerequisites:
  - $REFGENIE set to a refgenie config with hg38 and human_rDNA assets
  - All bioinformatics tools installed (cutadapt, fastp, seqtk, fastx,
    seqkit, fqdedup, bowtie2, samtools, bedtools, fastqc, preseq, R, etc.)
  - RUN_INTEGRATION_TESTS=true environment variable set

Run from the repository root:
    RUN_INTEGRATION_TESTS=true pytest tests/test_integration.py -v

To run a single scenario:
    RUN_INTEGRATION_TESTS=true pytest tests/test_integration.py -v -k se_basic

Keep output directories for debugging:
    KEEP_TEST_OUTPUTS=true RUN_INTEGRATION_TESTS=true pytest tests/test_integration.py -v
"""

import glob
import os
import shutil
import subprocess
import yaml
import pytest

# ---------------------------------------------------------------------------
# Gate: skip all integration tests unless explicitly enabled
# ---------------------------------------------------------------------------

INTEGRATION_ENABLED = os.environ.get("RUN_INTEGRATION_TESTS", "").lower() in (
    "1", "true", "yes"
)
KEEP_TEST_OUTPUTS = os.environ.get("KEEP_TEST_OUTPUTS", "").lower() in (
    "1", "true", "yes"
)

pytestmark = pytest.mark.skipif(
    not INTEGRATION_ENABLED,
    reason="Set RUN_INTEGRATION_TESTS=true to run integration tests",
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOOPER_CFG_DIR = os.path.join(REPO_ROOT, "tests", "looper_configs")

CORE_STATS = [
    "Raw_reads",
    "Trimmed_reads_R1",
    "Aligned_reads",
    "Alignment_rate",
    "Mapped_reads",
    "NRF",
    "PBC1",
    "PBC2",
]


def run_looper(looper_cfg, recover=False):
    """
    Run looper for a given looper config using the local compute package
    so the pipeline executes inline (not submitted to a job scheduler).
    Returns a CompletedProcess instance.
    """
    cmd = ["looper", "run", "-c", looper_cfg, "-p", "local"]
    if recover:
        # looper 2.1+ removed --recover; pass pypiper's -R flag via --command-extra.
        # Use --command-extra=-R (not -x -R) so argparse doesn't mistake -R for a flag.
        cmd.append("--command-extra=-R")
    return subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )


def assert_pipeline_succeeded(result, sample_dir):
    """Assert looper exited 0 and the pipeline log shows no failure."""
    assert result.returncode == 0, (
        f"looper exited {result.returncode}\n"
        f"STDOUT:\n{result.stdout[-3000:]}\n"
        f"STDERR:\n{result.stderr[-3000:]}"
    )
    log = os.path.join(sample_dir, "PEPPRO_log.md")
    if os.path.exists(log):
        assert "Pipeline failed" not in open(log).read(), (
            f"Pipeline log indicates failure:\n{open(log).read()[-2000:]}"
        )


def load_stats(sample_dir):
    """Load pipestat stats.yaml for a completed sample run.

    pipestat writes: PEPPRO.sample.<record_identifier>.<metric>
    Returns the flat metrics dict for the sample.
    """
    path = os.path.join(sample_dir, "stats.yaml")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        data = yaml.safe_load(f) or {}
    # Navigate PEPPRO -> sample -> <record_identifier>
    sample_name = os.path.basename(sample_dir)
    try:
        return data["PEPPRO"]["sample"][sample_name]
    except (KeyError, TypeError):
        return data


def assert_stats_keys(stats, expected_keys):
    """Assert that all expected metric keys are present in stats.yaml."""
    missing = [k for k in expected_keys if k not in stats]
    assert not missing, f"Missing stats keys: {missing}"


def scenario_output_dir(scenario):
    """Return the output_dir for a scenario (matches looper config's output_dir)."""
    return os.path.expandvars(f"$HOME/peppro_test_{scenario}/results_pipeline")


# ===========================================================================
# Base class: runs looper once per test class, shared by all test methods
# ===========================================================================

class PepproIntegrationBase:
    """
    Base class for single-scenario integration tests.

    setup_class runs looper once; all test methods in the class share the
    resulting output directory and CompletedProcess (self.result, self.sample_dir).
    teardown_class removes the output tree unless KEEP_TEST_OUTPUTS=true.
    """

    SCENARIO = None  # override in subclass, e.g. "se_basic"
    SAMPLE = None    # sample_name in the PEP CSV, e.g. "se_basic"

    @classmethod
    def setup_class(cls):
        cls.output_dir = scenario_output_dir(cls.SCENARIO)
        cls.sample_dir = os.path.join(cls.output_dir, cls.SAMPLE)
        # Create output_dir before running looper so it can write the
        # pipestat_config_PEPPRO.yaml file there before the pipeline starts.
        os.makedirs(cls.output_dir, exist_ok=True)
        cfg = os.path.join(LOOPER_CFG_DIR, f".looper_{cls.SCENARIO}.yaml")
        cls.result = run_looper(cfg)

    @classmethod
    def teardown_class(cls):
        if not KEEP_TEST_OUTPUTS:
            # Remove ${HOME}/peppro_test_<scenario}/ entirely
            parent = os.path.dirname(cls.output_dir)
            shutil.rmtree(parent, ignore_errors=True)


# ===========================================================================
# Integration test scenarios
# ===========================================================================

class Test_se_basic(PepproIntegrationBase):
    """SE PRO-seq baseline — cutadapt + seqtk, no UMI."""

    SCENARIO = "se_basic"
    SAMPLE   = "se_basic"

    def test_pipeline_runs(self):
        assert_pipeline_succeeded(self.result, self.sample_dir)

    def test_core_outputs_exist(self):
        assert glob.glob(os.path.join(self.sample_dir, "stats.yaml")), \
            f"stats.yaml missing in {self.sample_dir}"

    def test_core_stats_reported(self):
        assert_stats_keys(load_stats(self.sample_dir), CORE_STATS)

    def test_tss_score_reported(self):
        assert "TSS_coding_score" in load_stats(self.sample_dir)


class Test_pe_basic(PepproIntegrationBase):
    """PE PRO-seq baseline — cutadapt + seqtk, no UMI."""

    SCENARIO = "pe_basic"
    SAMPLE   = "pe_basic"

    def test_pipeline_runs(self):
        assert_pipeline_succeeded(self.result, self.sample_dir)

    def test_r2_trimmed_stats_reported(self):
        """PE-specific: Trimmed_reads_R2 and Trim_loss_rate_R2 must be present."""
        assert_stats_keys(load_stats(self.sample_dir), ["Trimmed_reads_R2", "Trim_loss_rate_R2"])

    def test_adapter_insertion_plot(self):
        """PE-only: adapter insertion distribution PDF should be generated."""
        pdfs = glob.glob(os.path.join(self.sample_dir, "cutadapt", "*insertion*distribution*.pdf"))
        assert pdfs, "Adapter insertion distribution plot not generated"

    def test_fastqc_r2_report(self):
        """PE-only: FastQC report for R2 should be generated."""
        r2_qc = glob.glob(os.path.join(self.sample_dir, "fastqc", "*R2*_fastqc.html"))
        assert r2_qc, "FastQC R2 report not generated"


class Test_se_groseq(PepproIntegrationBase):
    """SE GRO-seq protocol."""

    SCENARIO = "se_groseq"
    SAMPLE   = "se_groseq"

    def test_pipeline_runs(self):
        assert_pipeline_succeeded(self.result, self.sample_dir)

    def test_groseq_in_log(self):
        log = os.path.join(self.sample_dir, "PEPPRO_log.md")
        assert "GRO" in open(log).read()


class Test_se_umi(PepproIntegrationBase):
    """SE PRO-seq with 8-nt UMI, dedup with seqkit."""

    SCENARIO = "se_umi"
    SAMPLE   = "se_umi"

    def test_pipeline_runs(self):
        assert_pipeline_succeeded(self.result, self.sample_dir)

    def test_duplicate_reads_reported(self):
        assert "Duplicate_reads" in load_stats(self.sample_dir)


class Test_pe_umi(PepproIntegrationBase):
    """PE PRO-seq with 8-nt UMI, dedup with seqkit."""

    SCENARIO = "pe_umi"
    SAMPLE   = "pe_umi"

    def test_pipeline_runs(self):
        assert_pipeline_succeeded(self.result, self.sample_dir)

    def test_r2_trimmed_stats_reported(self):
        assert_stats_keys(load_stats(self.sample_dir), ["Trimmed_reads_R2", "Duplicate_reads"])


class Test_se_fastp(PepproIntegrationBase):
    """SE PRO-seq with fastp adapter removal."""

    SCENARIO = "se_fastp"
    SAMPLE   = "se_fastp"

    def test_pipeline_runs(self):
        assert_pipeline_succeeded(self.result, self.sample_dir)

    def test_fastp_report_generated(self):
        reports = glob.glob(os.path.join(self.sample_dir, "fastp", "*.html"))
        assert reports, "fastp HTML report not generated"


class Test_se_fastx(PepproIntegrationBase):
    """SE PRO-seq with fastx_trimmer."""

    SCENARIO = "se_fastx"
    SAMPLE   = "se_fastx"

    def test_pipeline_runs(self):
        assert_pipeline_succeeded(self.result, self.sample_dir)

    def test_core_stats_reported(self):
        assert_stats_keys(load_stats(self.sample_dir), CORE_STATS)


@pytest.mark.skipif(
    shutil.which("fqdedup") is None,
    reason="fqdedup not installed",
)
class Test_se_fqdedup(PepproIntegrationBase):
    """SE PRO-seq with 8-nt UMI, dedup with fqdedup."""

    SCENARIO = "se_fqdedup"
    SAMPLE   = "se_fqdedup"

    def test_pipeline_runs(self):
        assert_pipeline_succeeded(self.result, self.sample_dir)

    def test_duplicate_reads_reported(self):
        assert "Duplicate_reads" in load_stats(self.sample_dir)


class Test_se_scale(PepproIntegrationBase):
    """SE PRO-seq with --scale flag."""

    SCENARIO = "se_scale"
    SAMPLE   = "se_scale"

    def test_pipeline_runs(self):
        assert_pipeline_succeeded(self.result, self.sample_dir)

    def test_scale_stats_reported(self):
        """--scale passes through to bamSitesToWig; verify pipeline completes and reports stats."""
        assert_stats_keys(load_stats(self.sample_dir), CORE_STATS)


class Test_se_no_complexity(PepproIntegrationBase):
    """SE PRO-seq with --no-complexity (skip preseq)."""

    SCENARIO = "se_no_complexity"
    SAMPLE   = "se_no_complexity"

    def test_pipeline_runs(self):
        assert_pipeline_succeeded(self.result, self.sample_dir)

    def test_no_preseq_output(self):
        """Library complexity output should NOT be generated when skipped."""
        preseq = glob.glob(os.path.join(self.sample_dir, "QC_hg38", "*preseq*"))
        assert not preseq, "Preseq output unexpectedly generated with --no-complexity"


class Test_se_nofifo(PepproIntegrationBase):
    """SE PRO-seq with --noFIFO (disable named pipes in prealignments)."""

    SCENARIO = "se_nofifo"
    SAMPLE   = "se_nofifo"

    def test_pipeline_runs(self):
        assert_pipeline_succeeded(self.result, self.sample_dir)

    def test_core_stats_reported(self):
        assert_stats_keys(load_stats(self.sample_dir), CORE_STATS)


class Test_se_coverage(PepproIntegrationBase):
    """SE PRO-seq with --coverage flag."""

    SCENARIO = "se_coverage"
    SAMPLE   = "se_coverage"

    def test_pipeline_runs(self):
        assert_pipeline_succeeded(self.result, self.sample_dir)

    def test_core_stats_reported(self):
        assert_stats_keys(load_stats(self.sample_dir), CORE_STATS)


# ===========================================================================
# Recovery regression tests
# Uses a dedicated looper config (.looper_se_recovery.yaml) with its own
# output directory so it doesn't conflict with Test_se_basic.
# ===========================================================================

class Test_recovery:
    """
    Checkpoint recovery regression tests.

    Runs the pipeline twice against a dedicated output directory and verifies
    that completed pipeline stages are correctly skipped on recovery.
    """

    SCENARIO   = "se_recovery"
    SAMPLE     = "se_basic"   # reuses se_basic PEP config, same sample name
    OUTPUT_DIR = os.path.expandvars("$HOME/peppro_test_se_recovery")

    @classmethod
    def setup_class(cls):
        os.makedirs(cls.OUTPUT_DIR, exist_ok=True)
        cls.sample_dir = os.path.join(cls.OUTPUT_DIR, "results_pipeline", cls.SAMPLE)

    @classmethod
    def teardown_class(cls):
        if not KEEP_TEST_OUTPUTS:
            shutil.rmtree(cls.OUTPUT_DIR, ignore_errors=True)

    def test_recovery(self):
        """
        Full recovery regression:
        1. Run the pipeline to completion.
        2. Re-run with --recover; verify stages are skipped (log shows "Target exists").
        3. Delete processed_R1.flag to simulate a stale checkpoint; recovery should
           still succeed without the 'unmap_R1.fq not found' error that was fixed.
        """
        cfg = os.path.join(LOOPER_CFG_DIR, f".looper_{self.SCENARIO}.yaml")

        # --- Run 1: full pipeline ---
        r1 = run_looper(cfg)
        assert_pipeline_succeeded(r1, self.sample_dir)

        # --- Run 2: recover — verify checkpoint skipping ---
        r2 = run_looper(cfg, recover=True)
        assert_pipeline_succeeded(r2, self.sample_dir)
        log_content = open(os.path.join(self.sample_dir, "PEPPRO_log.md")).read()
        assert "Target exists" in log_content, \
            "Recovery run did not skip any completed stages"

        # --- Run 3: stale R1 flag regression ---
        # Delete processed_R1.flag so recovery must reconstruct unmap_fq1 path
        r1_flag = os.path.join(self.sample_dir, "fastq", "processed_R1.flag")
        if os.path.exists(r1_flag):
            os.remove(r1_flag)

        r3 = run_looper(cfg, recover=True)
        assert_pipeline_succeeded(r3, self.sample_dir)
        log_content = open(os.path.join(self.sample_dir, "PEPPRO_log.md")).read()
        assert "unmap_R1.fq" not in log_content, \
            "Recovery failed with missing unmap_R1.fq (regression)"
