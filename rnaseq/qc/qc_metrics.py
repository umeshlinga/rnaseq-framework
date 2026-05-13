"""
qc_metrics.py
-------------
QC Metrics Aggregation Module

Collects per-sample QC statistics across all pipeline steps
and generates a structured QC summary report.
"""

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Dict, Optional

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


@dataclass
class SampleQCSummary:
    """Aggregated QC metrics for a single sample across all pipeline steps."""
    sample_id: str

    # Trimming
    input_reads: int = 0
    surviving_reads: int = 0
    survival_rate: float = 0.0

    # Alignment
    uniquely_mapped: int = 0
    mapping_rate: float = 0.0
    multi_mapped: int = 0

    # Quantification
    assigned_fragments: int = 0
    assignment_rate: float = 0.0

    # Flags
    flagged: bool = False
    warnings: List[str] = field(default_factory=list)

    # QC thresholds
    MIN_SURVIVAL_RATE: float = 0.70
    MIN_MAPPING_RATE: float = 0.75
    MIN_ASSIGNMENT_RATE: float = 0.60

    def evaluate(self):
        """Check metrics against QC thresholds and set warnings."""
        if self.survival_rate < self.MIN_SURVIVAL_RATE:
            self.warnings.append(
                f"Low read survival rate: {self.survival_rate*100:.1f}% "
                f"(min: {self.MIN_SURVIVAL_RATE*100:.0f}%)"
            )
        if self.mapping_rate < self.MIN_MAPPING_RATE:
            self.warnings.append(
                f"Low mapping rate: {self.mapping_rate*100:.1f}% "
                f"(min: {self.MIN_MAPPING_RATE*100:.0f}%)"
            )
        if self.assignment_rate < self.MIN_ASSIGNMENT_RATE:
            self.warnings.append(
                f"Low feature assignment rate: {self.assignment_rate*100:.1f}% "
                f"(min: {self.MIN_ASSIGNMENT_RATE*100:.0f}%)"
            )
        self.flagged = len(self.warnings) > 0

    def to_dict(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "trimming": {
                "input_reads": self.input_reads,
                "surviving_reads": self.surviving_reads,
                "survival_rate_pct": round(self.survival_rate * 100, 2),
            },
            "alignment": {
                "uniquely_mapped": self.uniquely_mapped,
                "mapping_rate_pct": round(self.mapping_rate * 100, 2),
                "multi_mapped": self.multi_mapped,
            },
            "quantification": {
                "assigned_fragments": self.assigned_fragments,
                "assignment_rate_pct": round(self.assignment_rate * 100, 2),
            },
            "qc_status": "FAIL" if self.flagged else "PASS",
            "warnings": self.warnings,
        }


class QCMetricsAggregator:
    """
    Aggregates QC metrics across all pipeline steps for all samples.

    Collects statistics from trimming, alignment, and quantification results,
    applies QC thresholds, and generates a structured HTML + JSON report.

    Parameters
    ----------
    output_dir : str
        Directory to write QC reports.
    """

    def __init__(self, output_dir: str = "data/results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.samples: List[SampleQCSummary] = []

    def add_sample(
        self,
        sample_id: str,
        trimming_result=None,
        alignment_result=None,
        quant_result=None,
    ) -> SampleQCSummary:
        """
        Add QC metrics for a single sample.

        Parameters
        ----------
        sample_id : str
        trimming_result : TrimmingResult, optional
        alignment_result : AlignmentResult, optional
        quant_result : QuantificationResult, optional

        Returns
        -------
        SampleQCSummary
        """
        summary = SampleQCSummary(sample_id=sample_id)

        if trimming_result:
            summary.input_reads = trimming_result.input_reads
            summary.surviving_reads = trimming_result.surviving_reads
            summary.survival_rate = (
                trimming_result.surviving_reads / trimming_result.input_reads
                if trimming_result.input_reads > 0 else 0.0
            )

        if alignment_result:
            summary.uniquely_mapped = alignment_result.uniquely_mapped
            summary.mapping_rate = alignment_result.mapping_rate
            summary.multi_mapped = alignment_result.multi_mapped

        if quant_result:
            summary.assigned_fragments = quant_result.assigned_fragments
            summary.assignment_rate = quant_result.assignment_rate

        summary.evaluate()
        self.samples.append(summary)

        if summary.flagged:
            logger.warning(f"[{sample_id}] QC FAIL: {summary.warnings}")
        else:
            logger.info(f"[{sample_id}] QC PASS")

        return summary

    def write_json_report(self) -> str:
        """Write QC summary as a JSON file."""
        report_path = str(self.output_dir / "qc_summary.json")
        passed = sum(1 for s in self.samples if not s.flagged)
        report = {
            "total_samples": len(self.samples),
            "passed": passed,
            "failed": len(self.samples) - passed,
            "samples": [s.to_dict() for s in self.samples],
        }
        with open(report_path, "w") as f:
            json.dump(report, f, indent=2)
        logger.info(f"QC JSON report written: {report_path}")
        return report_path

    def write_html_report(self) -> str:
        """Generate a styled HTML QC report table."""
        report_path = str(self.output_dir / "qc_report.html")
        rows = ""
        for s in self.samples:
            status_color = "#d32f2f" if s.flagged else "#2e7d32"
            status_text = "FAIL" if s.flagged else "PASS"
            warnings = "<br>".join(s.warnings) if s.warnings else "—"
            rows += f"""
            <tr>
                <td>{s.sample_id}</td>
                <td>{s.input_reads:,}</td>
                <td>{s.surviving_reads:,}</td>
                <td>{s.survival_rate*100:.1f}%</td>
                <td>{s.uniquely_mapped:,}</td>
                <td>{s.mapping_rate*100:.1f}%</td>
                <td>{s.assigned_fragments:,}</td>
                <td>{s.assignment_rate*100:.1f}%</td>
                <td style="color:{status_color};font-weight:bold">{status_text}</td>
                <td style="font-size:12px;color:#666">{warnings}</td>
            </tr>"""

        html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>RNA-seq QC Report</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 40px; color: #333; }}
  h1 {{ font-size: 24px; margin-bottom: 4px; }}
  .summary {{ display: flex; gap: 24px; margin: 16px 0 24px; }}
  .card {{ background: #f5f5f5; padding: 16px 24px; border-radius: 8px; text-align: center; }}
  .card .num {{ font-size: 32px; font-weight: bold; }}
  .card .label {{ font-size: 13px; color: #666; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
  th {{ background: #1565c0; color: white; padding: 10px 12px; text-align: left; }}
  td {{ padding: 8px 12px; border-bottom: 1px solid #eee; }}
  tr:nth-child(even) {{ background: #fafafa; }}
</style>
</head>
<body>
<h1>RNA-seq QC Report</h1>
<p style="color:#666;font-size:13px">Generated by rnaseq-framework</p>
<div class="summary">
  <div class="card"><div class="num">{len(self.samples)}</div><div class="label">Total samples</div></div>
  <div class="card"><div class="num" style="color:#2e7d32">{sum(1 for s in self.samples if not s.flagged)}</div><div class="label">Passed QC</div></div>
  <div class="card"><div class="num" style="color:#d32f2f">{sum(1 for s in self.samples if s.flagged)}</div><div class="label">Failed QC</div></div>
</div>
<table>
<thead>
  <tr>
    <th>Sample</th>
    <th>Input reads</th>
    <th>Surviving reads</th>
    <th>Survival %</th>
    <th>Uniquely mapped</th>
    <th>Mapping %</th>
    <th>Assigned fragments</th>
    <th>Assignment %</th>
    <th>QC status</th>
    <th>Warnings</th>
  </tr>
</thead>
<tbody>{rows}</tbody>
</table>
</body>
</html>"""

        with open(report_path, "w") as f:
            f.write(html)
        logger.info(f"QC HTML report written: {report_path}")
        return report_path
