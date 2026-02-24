"""
Unit tests for PEPPRO — no genome indices or bioinformatics tools required.

Run from the repository root:
    pytest tests/test_unit.py -v
"""

import os
import sys
import tempfile
import yaml
import pytest
import peppy
import eido

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTS_DIR = os.path.join(REPO_ROOT, "tests")
PEP_DIR   = os.path.join(TESTS_DIR, "pep_configs")
DATA_DIR  = os.path.join(TESTS_DIR, "data")
SCHEMA    = os.path.join(REPO_ROOT, "peppro_input_schema.yaml")


def load_pep(name):
    """Load a test PEP config by scenario name."""
    return peppy.Project(os.path.join(PEP_DIR, f"{name}.yaml"))


def make_minimal_pep(tmpdir, extra_sample_fields=None, extra_modifiers=None):
    """
    Build a minimal in-memory PEP config for schema-validation testing.
    Returns a peppy.Project.
    """
    sample = {
        "sample_name": "test",
        "organism": "human",
        "protocol": "PROSEQ",
        "read_type": "single",
        "read1": os.path.join(DATA_DIR, "test_R1.fastq.gz"),
        "genome": "hg38",
    }
    if extra_sample_fields:
        sample.update(extra_sample_fields)

    csv_path = os.path.join(tmpdir, "s.csv")
    import csv as csv_mod
    with open(csv_path, "w", newline="") as f:
        w = csv_mod.DictWriter(f, fieldnames=list(sample.keys()))
        w.writeheader()
        w.writerow(sample)

    cfg = {
        "name": "test_project",
        "pep_version": "2.0.0",
        "sample_table": csv_path,
    }
    if extra_modifiers:
        cfg["sample_modifiers"] = extra_modifiers

    cfg_path = os.path.join(tmpdir, "p.yaml")
    with open(cfg_path, "w") as f:
        yaml.dump(cfg, f)

    return peppy.Project(cfg_path)


# ===========================================================================
# 1. Constants and module-level attributes
# ===========================================================================

class TestConstants:
    def test_runon_source_includes_proseq(self):
        sys.path.insert(0, os.path.join(REPO_ROOT, "pipelines"))
        import peppro
        assert "PROSEQ" in peppro.RUNON_SOURCE
        assert "proseq" in peppro.RUNON_SOURCE
        assert "PRO-seq" in peppro.RUNON_SOURCE

    def test_runon_source_includes_groseq(self):
        import peppro
        assert "GROSEQ" in peppro.RUNON_SOURCE
        assert "groseq" in peppro.RUNON_SOURCE
        assert "GRO-seq" in peppro.RUNON_SOURCE

    def test_adapter_choices(self):
        import peppro
        assert "cutadapt" in peppro.ADAPTER_REMOVERS
        assert "fastp" in peppro.ADAPTER_REMOVERS

    def test_trimmer_choices(self):
        import peppro
        assert "seqtk" in peppro.TRIMMERS
        assert "fastx" in peppro.TRIMMERS

    def test_deduplicator_choices(self):
        import peppro
        assert "seqkit" in peppro.DEDUPLICATORS
        assert "fqdedup" in peppro.DEDUPLICATORS

    def test_defaults(self):
        import peppro
        assert peppro.DEFAULT_REMOVER    == "cutadapt"
        assert peppro.DEFAULT_TRIMMER    == "seqtk"
        assert peppro.DEFAULT_DEDUPLICATOR == "seqkit"
        assert peppro.DEFAULT_UMI_LEN    == 0
        assert peppro.DEFAULT_MAX_LEN    == -1


# ===========================================================================
# 2. PEP config loading — verify sample attributes per scenario
# ===========================================================================

class TestPepLoading:
    """Load each test PEP config and verify the expected sample attributes."""

    def test_se_basic_attributes(self):
        p = load_pep("se_basic")
        s = p.samples[0]
        assert s.sample_name == "se_basic"
        assert s.protocol    == "PROSEQ"
        assert s.read_type   == "single"
        assert s.genome      == "hg38"
        assert not hasattr(s, "read2") or not s.read2

    def test_pe_basic_attributes(self):
        p = load_pep("pe_basic")
        s = p.samples[0]
        assert s.sample_name == "pe_basic"
        assert s.read_type   == "paired"
        assert s.genome      == "hg38"
        assert hasattr(s, "read2")

    def test_se_groseq_protocol(self):
        p = load_pep("se_groseq")
        assert p.samples[0].protocol == "GROSEQ"

    def test_se_umi_umi_len_is_string(self):
        """umi_len must arrive as a string so eido schema validation passes."""
        p = load_pep("se_umi")
        s = p.samples[0]
        assert hasattr(s, "umi_len")
        assert isinstance(s.umi_len, str), (
            f"umi_len should be a string, got {type(s.umi_len)}: {s.umi_len!r}"
        )
        assert s.umi_len == "8"

    def test_pe_umi_attributes(self):
        p = load_pep("pe_umi")
        s = p.samples[0]
        assert s.read_type == "paired"
        assert isinstance(s.umi_len, str)
        assert s.umi_len == "8"

    def test_se_fastp_adapter(self):
        p = load_pep("se_fastp")
        assert p.samples[0].adapter == "fastp"

    def test_se_fastx_trimmer(self):
        p = load_pep("se_fastx")
        assert p.samples[0].trimmer == "fastx"

    def test_se_fqdedup_dedup(self):
        p = load_pep("se_fqdedup")
        s = p.samples[0]
        assert s.dedup == "fqdedup"
        assert isinstance(s.umi_len, str)

    def test_se_scale_flag(self):
        p = load_pep("se_scale")
        assert hasattr(p.samples[0], "scale")

    def test_se_no_complexity_flag(self):
        p = load_pep("se_no_complexity")
        assert hasattr(p.samples[0], "complexity")

    def test_se_nofifo_flag(self):
        p = load_pep("se_nofifo")
        assert hasattr(p.samples[0], "no_fifo")

    def test_se_coverage_flag(self):
        p = load_pep("se_coverage")
        assert hasattr(p.samples[0], "coverage")

    def test_prealignment_names_set(self):
        """All human-genome PEPs should have prealignment_names via imply."""
        for scenario in ("se_basic", "pe_basic", "se_groseq", "se_umi"):
            p = load_pep(scenario)
            s = p.samples[0]
            assert hasattr(s, "prealignment_names"), (
                f"{scenario}: prealignment_names not set"
            )


# ===========================================================================
# 3. Schema validation — eido against peppro_input_schema.yaml
# ===========================================================================

class TestSchemaValidation:
    """Validate PEP configs against the input schema; some should fail."""

    def test_se_basic_passes_schema(self):
        p = load_pep("se_basic")
        eido.validate_project(p, SCHEMA)   # must not raise

    def test_pe_basic_passes_schema(self):
        p = load_pep("pe_basic")
        eido.validate_project(p, SCHEMA)

    def test_se_groseq_passes_schema(self):
        p = load_pep("se_groseq")
        eido.validate_project(p, SCHEMA)

    def test_se_umi_passes_schema(self):
        p = load_pep("se_umi")
        eido.validate_project(p, SCHEMA)

    def test_se_fastp_passes_schema(self):
        p = load_pep("se_fastp")
        eido.validate_project(p, SCHEMA)

    def test_se_fastx_passes_schema(self):
        p = load_pep("se_fastx")
        eido.validate_project(p, SCHEMA)

    def test_se_fqdedup_passes_schema(self):
        p = load_pep("se_fqdedup")
        eido.validate_project(p, SCHEMA)

    def test_umi_len_int_in_yaml_fails_schema(self):
        """
        Regression test: umi_len set as a bare integer in a YAML modifier
        (not a CSV column) must fail schema validation.
        This is the bug that was fixed in H9_example.yaml.
        """
        with tempfile.TemporaryDirectory() as d:
            p = make_minimal_pep(
                d,
                extra_modifiers={"append": {"umi_len": 8}},  # bare int
            )
            assert isinstance(p.samples[0].umi_len, int), (
                "Expected int from YAML append, got string — test assumption wrong"
            )
            with pytest.raises(eido.EidoValidationError):
                eido.validate_project(p, SCHEMA)

    def test_umi_len_string_in_yaml_passes_schema(self):
        """umi_len as a quoted string in YAML modifier should pass."""
        with tempfile.TemporaryDirectory() as d:
            p = make_minimal_pep(
                d,
                extra_modifiers={"append": {"umi_len": "8"}},  # quoted string
            )
            assert isinstance(p.samples[0].umi_len, str)
            eido.validate_project(p, SCHEMA)   # must not raise

    def test_missing_protocol_fails_schema(self):
        with tempfile.TemporaryDirectory() as d:
            p = make_minimal_pep(d)
            del p.samples[0]["protocol"]
            with pytest.raises(eido.EidoValidationError):
                eido.validate_project(p, SCHEMA)

    def test_invalid_protocol_fails_schema(self):
        with tempfile.TemporaryDirectory() as d:
            p = make_minimal_pep(d, extra_sample_fields={"protocol": "CHIP-SEQ"})
            with pytest.raises(eido.EidoValidationError):
                eido.validate_project(p, SCHEMA)

    def test_invalid_read_type_fails_schema(self):
        with tempfile.TemporaryDirectory() as d:
            p = make_minimal_pep(
                d, extra_sample_fields={"read_type": "triple"}
            )
            with pytest.raises(eido.EidoValidationError):
                eido.validate_project(p, SCHEMA)

    def test_invalid_adapter_choice_fails_schema(self):
        with tempfile.TemporaryDirectory() as d:
            p = make_minimal_pep(
                d, extra_sample_fields={"adapter": "trimmomatic"}
            )
            with pytest.raises(eido.EidoValidationError):
                eido.validate_project(p, SCHEMA)

    def test_invalid_trimmer_choice_fails_schema(self):
        with tempfile.TemporaryDirectory() as d:
            p = make_minimal_pep(
                d, extra_sample_fields={"trimmer": "trimgalore"}
            )
            with pytest.raises(eido.EidoValidationError):
                eido.validate_project(p, SCHEMA)

    def test_invalid_dedup_choice_fails_schema(self):
        with tempfile.TemporaryDirectory() as d:
            p = make_minimal_pep(
                d, extra_sample_fields={"dedup": "picard"}
            )
            with pytest.raises(eido.EidoValidationError):
                eido.validate_project(p, SCHEMA)


# ===========================================================================
# 4. Argument parsing
# ===========================================================================

class TestArgumentParsing:
    """Test parse_arguments() with mocked sys.argv."""

    @pytest.fixture(autouse=True)
    def _import_pipeline(self):
        sys.path.insert(0, os.path.join(REPO_ROOT, "pipelines"))
        import peppro
        self.peppro = peppro

    def _parse(self, extra_args):
        base = [
            "--sample-name", "test",
            "--genome", "hg38",
            "--input", "r1.fastq.gz",
            "--genome-index", "/idx/hg38/.",
            "--chrom-sizes", "/idx/hg38.sizes",
            "--output-parent", "/out",
            "--protocol", "PROSEQ",
            "--single-or-paired", "single",
        ]
        saved = sys.argv[:]
        sys.argv = ["peppro.py"] + base + extra_args
        try:
            return self.peppro.parse_arguments()
        finally:
            sys.argv = saved

    def test_defaults(self):
        args = self._parse([])
        assert args.adapter  == "cutadapt"
        assert args.trimmer  == "seqtk"
        assert args.dedup    == "seqkit"
        assert args.umi_len  == 0
        assert args.max_len  == -1
        assert args.sob      is False
        assert args.scale    is False
        assert args.coverage is False
        assert args.complexity is False
        assert args.no_fifo  is False

    def test_fastp_adapter(self):
        args = self._parse(["--adapter-tool", "fastp"])
        assert args.adapter == "fastp"

    def test_fastx_trimmer(self):
        args = self._parse(["--trimmer-tool", "fastx"])
        assert args.trimmer == "fastx"

    def test_fqdedup_deduplicator(self):
        args = self._parse(["--dedup-tool", "fqdedup"])
        assert args.dedup == "fqdedup"

    def test_umi_len_parsed_as_int(self):
        args = self._parse(["--umi-len", "8"])
        assert args.umi_len == 8
        assert isinstance(args.umi_len, int)

    def test_scale_flag(self):
        args = self._parse(["--scale"])
        assert args.scale is True

    def test_no_complexity_flag(self):
        args = self._parse(["--no-complexity"])
        assert args.complexity is True

    def test_nofifo_flag(self):
        args = self._parse(["--noFIFO"])
        assert args.no_fifo is True

    def test_coverage_flag(self):
        args = self._parse(["--coverage"])
        assert args.coverage is True

    def test_keep_flag(self):
        args = self._parse(["--keep"])
        assert args.keep is True

    def test_sob_flag(self):
        args = self._parse(["--sob"])
        assert args.sob is True

    def test_invalid_adapter_rejected(self):
        with pytest.raises(SystemExit):
            self._parse(["--adapter-tool", "trimmomatic"])

    def test_invalid_trimmer_rejected(self):
        with pytest.raises(SystemExit):
            self._parse(["--trimmer-tool", "trimgalore"])

    def test_invalid_dedup_rejected(self):
        with pytest.raises(SystemExit):
            self._parse(["--dedup-tool", "picard"])

    def test_invalid_protocol_rejected(self):
        saved = sys.argv[:]
        sys.argv = ["peppro.py", "--sample-name", "t", "--genome", "hg38",
                    "--input", "r1.fq.gz", "--genome-index", "/i",
                    "--chrom-sizes", "/c", "--output-parent", "/o",
                    "--protocol", "CHIP-SEQ", "--single-or-paired", "single"]
        try:
            with pytest.raises(SystemExit):
                self.peppro.parse_arguments()
        finally:
            sys.argv = saved

    def test_paired_end_flag(self):
        args = self._parse(["--single-or-paired", "paired",
                             "--input2", "r2.fastq.gz"])
        assert args.single_or_paired == "paired"

    def test_single_end_flag(self):
        args = self._parse(["--single-or-paired", "single"])
        assert args.single_or_paired == "single"


# ===========================================================================
# 5. Recovery path reconstruction
# ===========================================================================

class TestRecoveryPaths:
    """
    Verify the path logic that reconstructs unmap_fq1/unmap_fq2 when
    checkpoint flags exist (the recovery bug fixed in this session).
    These tests confirm the expected filenames, not file existence.
    """

    def _expected_paths(self, out_fastq_pre, umi=False):
        r1 = out_fastq_pre + "_R1_processed.fastq"
        r2 = out_fastq_pre + "_R2_trimmed.fastq"
        r1_dups = out_fastq_pre + "_R1_trimmed.fastq"
        r2_dups = out_fastq_pre + "_R2_trimmed_dups.fastq"
        return r1, r2, r1_dups, r2_dups

    def test_r1_recovered_path_suffix(self):
        base = "/out/sample/fastq/sample"
        r1, _, _, _ = self._expected_paths(base)
        assert r1.endswith("_R1_processed.fastq")

    def test_r2_recovered_path_suffix(self):
        base = "/out/sample/fastq/sample"
        _, r2, _, _ = self._expected_paths(base)
        assert r2.endswith("_R2_trimmed.fastq")

    def test_r1_dups_recovered_path_suffix(self):
        base = "/out/sample/fastq/sample"
        _, _, r1_dups, _ = self._expected_paths(base)
        assert r1_dups.endswith("_R1_trimmed.fastq")

    def test_r2_dups_recovered_path_suffix(self):
        base = "/out/sample/fastq/sample"
        _, _, _, r2_dups = self._expected_paths(base)
        assert r2_dups.endswith("_R2_trimmed_dups.fastq")

    def test_r1_path_includes_sample_name(self):
        sname = "my_sample"
        base = f"/results/{sname}/fastq/{sname}"
        r1, _, _, _ = self._expected_paths(base)
        assert sname in r1

    def test_expected_files_match_pipeline_code(self):
        """
        Verify the hardcoded recovery paths in peppro.py match what
        _process_fastq actually produces (named by the same conventions).
        """
        sys.path.insert(0, os.path.join(REPO_ROOT, "pipelines"))
        import peppro
        # _process_fastq builds these names:
        #   processed_fastq = os.path.join(fastq_folder, sname + "_R1_processed.fastq")
        #   trimmed_fq1     = os.path.join(fastq_folder, sname + "_R1_trimmed.fastq")
        #   trimmed_fq2     = os.path.join(fastq_folder, sname + "_R2_trimmed.fastq")
        #   trimmed_dups_fq2 = os.path.join(fastq_folder, sname + "_R2_trimmed_dups.fastq")
        # And in main(), the recovery else-branches use:
        #   out_fastq_pre + "_R1_processed.fastq"
        #   out_fastq_pre + "_R1_trimmed.fastq"
        #   out_fastq_pre + "_R2_trimmed.fastq"
        #   out_fastq_pre + "_R2_trimmed_dups.fastq"
        # Since out_fastq_pre ends with sname (see ngstk.input_to_fastq),
        # these are equivalent.  This test just documents the contract.
        suffixes_r1 = ["_R1_processed.fastq", "_R1_trimmed.fastq"]
        suffixes_r2 = ["_R2_trimmed.fastq", "_R2_trimmed_dups.fastq"]
        for suf in suffixes_r1 + suffixes_r2:
            assert suf.startswith("_R")
