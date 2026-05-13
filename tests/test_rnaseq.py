"""
test_rnaseq.py
--------------
Unit tests for the RNA-seq framework modules.
"""

import os
import json
import pytest
import tempfile
from unittest.mock import patch, MagicMock, mock_open
from pathlib import Path

from rnaseq.trimming.trimmer import ReadTrimmer, TrimmingResult
from rnaseq.quantification.quantifier import GeneQuantifier
from rnaseq.qc.qc_metrics import QCMetricsAggregator, SampleQCSummary


# ── Fixtures ──────────────────────────────────────────────────────

def make_counts_file(path: str, genes: dict):
    """Write a minimal featureCounts output file."""
    with open(path, "w") as f:
        f.write("# Program: featureCounts\n")
        f.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\tSample\n")
        for gene, count in genes.items():
            f.write(f"{gene}\tchr1\t1000\t2000\t+\t1000\t{count}\n")


def make_summary_file(path: str, assigned: int, multimapping: int, nofeatures: int):
    """Write a featureCounts summary file."""
    with open(path, "w") as f:
        f.write("Status\tSample\n")
        f.write(f"Assigned\t{assigned}\n")
        f.write(f"Unassigned_MultiMapping\t{multimapping}\n")
        f.write(f"Unassigned_NoFeatures\t{nofeatures}\n")


# ── Trimmer Tests ─────────────────────────────────────────────────

class TestReadTrimmer:

    def test_survival_rate_calculation(self):
        """TrimmingResult.survival_rate should compute correctly."""
        result = TrimmingResult(
            sample_id="s1",
            trimmed_r1="out.fastq.gz",
            input_reads=1000,
            surviving_reads=920,
            success=True,
        )
        assert result.survival_rate == 92.0

    def test_zero_input_survival_rate(self):
        """Zero input reads should give 0.0 survival rate (no division error)."""
        result = TrimmingResult(sample_id="s1", trimmed_r1="out.fastq.gz")
        assert result.survival_rate == 0.0

    def test_summary_keys(self):
        """TrimmingResult.summary should include all expected keys."""
        result = TrimmingResult(
            sample_id="s1",
            trimmed_r1="out.fastq.gz",
            input_reads=500,
            surviving_reads=480,
            success=True,
        )
        s = result.summary
        assert "sample_id" in s
        assert "survival_rate_pct" in s
        assert s["survival_rate_pct"] == 96.0

    @patch("rnaseq.trimming.trimmer.subprocess.run")
    def test_trim_paired_success(self, mock_run, tmp_path):
        """Paired trimming should parse Trimmomatic output correctly."""
        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stderr = (
            "TrimmomaticPE: Started\n"
            "Input Read Pairs: 1000 Both Surviving: 920 (92.00%) "
            "Forward Only Surviving: 40 (4.00%) "
            "Reverse Only Surviving: 20 (2.00%) "
            "Dropped: 20 (2.00%)\n"
        )
        mock_run.return_value = mock_proc

        trimmer = ReadTrimmer(
            adapters="adapters.fa",
            output_dir=str(tmp_path),
        )
        result = trimmer.trim_paired("R1.fastq.gz", "R2.fastq.gz", "sample1")

        assert result.success
        assert result.input_reads == 1000
        assert result.surviving_reads == 920
        assert result.dropped_reads == 20

    @patch("rnaseq.trimming.trimmer.subprocess.run")
    def test_trim_failure_handling(self, mock_run, tmp_path):
        """Failed Trimmomatic run should set success=False and error_message."""
        mock_proc = MagicMock()
        mock_proc.returncode = 1
        mock_proc.stderr = "Error: adapter file not found"
        mock_run.return_value = mock_proc

        trimmer = ReadTrimmer(adapters="missing.fa", output_dir=str(tmp_path))
        result = trimmer.trim_single("input.fastq.gz", "sample_fail")

        assert not result.success
        assert "Error" in result.error_message


# ── Quantifier Tests ──────────────────────────────────────────────

class TestGeneQuantifier:

    def test_build_count_matrix(self, tmp_path):
        """Count matrix should merge per-sample count files correctly."""
        counts1 = str(tmp_path / "s1.counts.txt")
        counts2 = str(tmp_path / "s2.counts.txt")
        make_counts_file(counts1, {"BRCA1": 100, "TP53": 250, "EGFR": 50})
        make_counts_file(counts2, {"BRCA1": 80, "TP53": 300, "EGFR": 120})

        quantifier = GeneQuantifier(gtf="dummy.gtf", output_dir=str(tmp_path))
        matrix = quantifier.build_count_matrix(
            [counts1, counts2], ["s1", "s2"]
        )

        assert "BRCA1" in matrix
        assert "TP53" in matrix
        assert matrix["BRCA1"]["s1"] == 100
        assert matrix["TP53"]["s2"] == 300

    def test_count_matrix_written_to_file(self, tmp_path):
        """build_count_matrix with output_path should write a TSV file."""
        counts1 = str(tmp_path / "s1.counts.txt")
        make_counts_file(counts1, {"GENE_A": 500, "GENE_B": 200})

        quantifier = GeneQuantifier(gtf="dummy.gtf", output_dir=str(tmp_path))
        out_path = str(tmp_path / "matrix.tsv")
        quantifier.build_count_matrix([counts1], ["sample1"], output_path=out_path)

        assert Path(out_path).exists()
        with open(out_path) as f:
            lines = f.readlines()
        assert lines[0].strip() == "gene_id\tsample1"
        assert len(lines) == 3  # header + 2 genes

    def test_rpkm_computation(self):
        """RPKM should normalise by gene length and library size."""
        matrix = {
            "GENE_A": {"s1": 1000},
            "GENE_B": {"s1": 500},
        }
        gene_lengths = {"GENE_A": 2000, "GENE_B": 1000}

        quantifier = GeneQuantifier(gtf="dummy.gtf")
        rpkm = quantifier.compute_rpkm(matrix, gene_lengths)

        assert "GENE_A" in rpkm
        assert "GENE_B" in rpkm
        # GENE_A: 1000 / (2.0kb * 1.5M reads) = 333.33
        # GENE_B: 500 / (1.0kb * 1.5M reads) = 333.33
        assert rpkm["GENE_A"]["s1"] > 0
        assert rpkm["GENE_B"]["s1"] > 0

    def test_rpkm_zero_length_fallback(self):
        """Genes missing from gene_lengths should use 1000bp default."""
        matrix = {"UNKNOWN_GENE": {"s1": 100}}
        quantifier = GeneQuantifier(gtf="dummy.gtf")
        rpkm = quantifier.compute_rpkm(matrix, gene_lengths={})
        assert rpkm["UNKNOWN_GENE"]["s1"] >= 0


# ── QC Metrics Tests ──────────────────────────────────────────────

class TestQCMetricsAggregator:

    def test_sample_passes_qc(self, tmp_path):
        """A sample with good metrics should pass QC with no warnings."""
        summary = SampleQCSummary(
            sample_id="good_sample",
            input_reads=10_000_000,
            surviving_reads=9_500_000,
            survival_rate=0.95,
            uniquely_mapped=8_800_000,
            mapping_rate=0.88,
            assigned_fragments=8_000_000,
            assignment_rate=0.84,
        )
        summary.evaluate()
        assert not summary.flagged
        assert len(summary.warnings) == 0

    def test_sample_fails_low_mapping(self):
        """Low mapping rate should trigger a QC warning and flag the sample."""
        summary = SampleQCSummary(
            sample_id="low_mapping",
            survival_rate=0.92,
            mapping_rate=0.50,  # below 0.75 threshold
            assignment_rate=0.75,
        )
        summary.evaluate()
        assert summary.flagged
        assert any("mapping" in w.lower() for w in summary.warnings)

    def test_sample_fails_low_survival(self):
        """Low read survival should trigger a QC warning."""
        summary = SampleQCSummary(
            sample_id="low_survival",
            survival_rate=0.55,  # below 0.70 threshold
            mapping_rate=0.85,
            assignment_rate=0.80,
        )
        summary.evaluate()
        assert summary.flagged
        assert any("survival" in w.lower() for w in summary.warnings)

    def test_write_json_report(self, tmp_path):
        """JSON report should contain expected structure."""
        agg = QCMetricsAggregator(output_dir=str(tmp_path))

        trim = MagicMock()
        trim.input_reads = 1_000_000
        trim.surviving_reads = 950_000

        align = MagicMock()
        align.uniquely_mapped = 880_000
        align.mapping_rate = 0.88
        align.multi_mapped = 20_000

        quant = MagicMock()
        quant.assigned_fragments = 820_000
        quant.assignment_rate = 0.82

        agg.add_sample("sample_A", trim, align, quant)
        report_path = agg.write_json_report()

        assert Path(report_path).exists()
        with open(report_path) as f:
            report = json.load(f)

        assert report["total_samples"] == 1
        assert report["passed"] == 1
        assert report["samples"][0]["sample_id"] == "sample_A"
        assert report["samples"][0]["qc_status"] == "PASS"

    def test_write_html_report(self, tmp_path):
        """HTML report should be valid HTML with expected content."""
        agg = QCMetricsAggregator(output_dir=str(tmp_path))

        trim = MagicMock()
        trim.input_reads = 500_000
        trim.surviving_reads = 480_000

        align = MagicMock()
        align.uniquely_mapped = 430_000
        align.mapping_rate = 0.80
        align.multi_mapped = 10_000

        quant = MagicMock()
        quant.assigned_fragments = 400_000
        quant.assignment_rate = 0.78

        agg.add_sample("sample_B", trim, align, quant)
        html_path = agg.write_html_report()

        assert Path(html_path).exists()
        with open(html_path) as f:
            content = f.read()
        assert "sample_B" in content
        assert "PASS" in content
        assert "<table>" in content

    def test_multiple_samples_mixed_results(self, tmp_path):
        """QC aggregator should correctly count pass/fail across samples."""
        agg = QCMetricsAggregator(output_dir=str(tmp_path))

        for i, (mapping_rate, assignment_rate) in enumerate([
            (0.90, 0.85),   # PASS
            (0.45, 0.80),   # FAIL - low mapping
            (0.85, 0.40),   # FAIL - low assignment
        ]):
            trim = MagicMock()
            trim.input_reads = 1_000_000
            trim.surviving_reads = 950_000
            align = MagicMock()
            align.uniquely_mapped = int(mapping_rate * 950_000)
            align.mapping_rate = mapping_rate
            align.multi_mapped = 10_000
            quant = MagicMock()
            quant.assigned_fragments = int(assignment_rate * 900_000)
            quant.assignment_rate = assignment_rate
            agg.add_sample(f"sample_{i}", trim, align, quant)

        report_path = agg.write_json_report()
        with open(report_path) as f:
            report = json.load(f)

        assert report["total_samples"] == 3
        assert report["passed"] == 1
        assert report["failed"] == 2
