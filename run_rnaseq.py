"""
run_rnaseq.py
-------------
RNA-seq Pipeline Orchestrator

End-to-end pipeline: Trimming → Alignment → Quantification → QC Report.
Reads a sample sheet CSV and processes each sample.

Usage
-----
    python run_rnaseq.py --config config/config.yaml --samples sample_sheet.csv
    python run_rnaseq.py --config config/config.yaml --samples sample_sheet.csv --threads 8
"""

import argparse
import csv
import json
import logging
import sys
import time
import yaml
from datetime import datetime
from pathlib import Path
from typing import List, Dict

from rnaseq.trimming.trimmer import ReadTrimmer
from rnaseq.alignment.rna_aligner import RNASeqAligner
from rnaseq.quantification.quantifier import GeneQuantifier
from rnaseq.qc.qc_metrics import QCMetricsAggregator

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(
            f"rnaseq_run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        ),
    ],
)
logger = logging.getLogger(__name__)


def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_sample_sheet(path: str) -> List[Dict]:
    samples = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            samples.append(row)
    logger.info(f"Loaded {len(samples)} samples from {path}")
    return samples


def run_sample(sample: Dict, config: dict, qc: QCMetricsAggregator) -> Dict:
    """
    Run the full RNA-seq pipeline for a single sample.

    Parameters
    ----------
    sample : dict
        Sample metadata: sample_id, fastq_r1, fastq_r2 (optional).
    config : dict
        Pipeline configuration.
    qc : QCMetricsAggregator
        Shared QC aggregator.

    Returns
    -------
    dict
        Per-sample results summary.
    """
    sample_id = sample["sample_id"]
    fastq_r1 = sample["fastq_r1"]
    fastq_r2 = sample.get("fastq_r2")
    paired = fastq_r2 is not None

    logger.info(f"{'='*60}")
    logger.info(f"Processing sample: {sample_id} ({'paired-end' if paired else 'single-end'})")
    logger.info(f"{'='*60}")

    start = time.time()
    results = {"sample_id": sample_id, "status": "failed", "steps": {}}

    trim_result = align_result = quant_result = None

    # ── Step 1: Adapter Trimming ──────────────────────────────────
    logger.info(f"[{sample_id}] Step 1/3: Read Trimming")
    try:
        trimmer = ReadTrimmer(
            adapters=config["trimming"]["adapters"],
            output_dir=config["output"]["trimmed_dir"],
            leading=config["trimming"]["leading"],
            trailing=config["trimming"]["trailing"],
            sliding_window=config["trimming"]["sliding_window"],
            min_length=config["trimming"]["min_length"],
            threads=config["alignment"]["threads"],
        )
        if paired:
            trim_result = trimmer.trim_paired(fastq_r1, fastq_r2, sample_id)
        else:
            trim_result = trimmer.trim_single(fastq_r1, sample_id)

        results["steps"]["trimming"] = trim_result.summary

        if not trim_result.success:
            logger.error(f"[{sample_id}] Trimming failed")
            return results

        if trim_result.survival_rate < config["trimming"]["min_survival_rate"] * 100:
            logger.warning(
                f"[{sample_id}] Low survival rate: {trim_result.survival_rate}%"
            )

    except Exception as e:
        logger.error(f"[{sample_id}] Trimming error: {e}")
        results["steps"]["trimming"] = {"error": str(e)}
        return results

    # ── Step 2: Alignment ─────────────────────────────────────────
    logger.info(f"[{sample_id}] Step 2/3: STAR Alignment")
    try:
        aligner = RNASeqAligner(
            genome_dir=config["reference"]["star_index"],
            output_dir=config["output"]["aligned_dir"],
            threads=config["alignment"]["threads"],
            two_pass=config["alignment"]["two_pass"],
        )
        align_result = aligner.align(
            fastq_r1=trim_result.trimmed_r1,
            sample_id=sample_id,
            fastq_r2=trim_result.trimmed_r2,
            gtf=config["reference"]["gtf"],
        )
        results["steps"]["alignment"] = align_result.summary

        if not align_result.success:
            logger.error(f"[{sample_id}] Alignment failed")
            return results

    except Exception as e:
        logger.error(f"[{sample_id}] Alignment error: {e}")
        results["steps"]["alignment"] = {"error": str(e)}
        return results

    # ── Step 3: Quantification ────────────────────────────────────
    logger.info(f"[{sample_id}] Step 3/3: Gene Quantification")
    try:
        quantifier = GeneQuantifier(
            gtf=config["reference"]["gtf"],
            output_dir=config["output"]["counts_dir"],
            paired_end=paired,
            strand_specific=config["quantification"]["strand_specific"],
            threads=config["alignment"]["threads"],
        )
        quant_result = quantifier.quantify(
            bam_path=align_result.bam_path,
            sample_id=sample_id,
        )
        results["steps"]["quantification"] = quant_result.summary

        if not quant_result.success:
            logger.error(f"[{sample_id}] Quantification failed")

    except Exception as e:
        logger.error(f"[{sample_id}] Quantification error: {e}")
        results["steps"]["quantification"] = {"error": str(e)}

    # ── Aggregate QC ──────────────────────────────────────────────
    qc.add_sample(
        sample_id=sample_id,
        trimming_result=trim_result,
        alignment_result=align_result,
        quant_result=quant_result,
    )

    elapsed = round(time.time() - start, 2)
    results["status"] = "completed"
    results["elapsed_seconds"] = elapsed
    logger.info(f"[{sample_id}] Completed in {elapsed}s")
    return results


def build_count_matrix(samples: List[Dict], config: dict) -> str:
    """Build the final count matrix from all per-sample count files."""
    counts_dir = Path(config["output"]["counts_dir"])
    sample_ids = [s["sample_id"] for s in samples]
    count_files = [str(counts_dir / f"{sid}.counts.txt") for sid in sample_ids]

    valid = [(f, sid) for f, sid in zip(count_files, sample_ids) if Path(f).exists()]
    if not valid:
        logger.warning("No count files found — skipping matrix build")
        return ""

    files, ids = zip(*valid)
    quantifier = GeneQuantifier(gtf=config["reference"]["gtf"])
    matrix_path = str(Path(config["output"]["results_dir"]) / "count_matrix.tsv")
    Path(config["output"]["results_dir"]).mkdir(parents=True, exist_ok=True)
    quantifier.build_count_matrix(list(files), list(ids), output_path=matrix_path)
    return matrix_path


def main():
    parser = argparse.ArgumentParser(description="RNA-seq Analysis Pipeline")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--samples", required=True, help="Path to sample sheet CSV")
    parser.add_argument("--threads", type=int, default=None, help="Override thread count")
    args = parser.parse_args()

    config = load_config(args.config)
    if args.threads:
        config["alignment"]["threads"] = args.threads

    samples = load_sample_sheet(args.samples)
    qc = QCMetricsAggregator(output_dir=config["output"]["results_dir"])
    all_results = []

    logger.info(f"Starting RNA-seq pipeline — {len(samples)} sample(s)")
    pipeline_start = time.time()

    for sample in samples:
        result = run_sample(sample, config, qc)
        all_results.append(result)

    # Build count matrix
    logger.info("Building count matrix...")
    matrix_path = build_count_matrix(samples, config)

    # Write QC reports
    json_report = qc.write_json_report()
    html_report = qc.write_html_report()

    # Write pipeline summary
    summary_path = Path(config["output"]["results_dir"]) / "pipeline_summary.json"
    completed = sum(1 for r in all_results if r["status"] == "completed")
    with open(summary_path, "w") as f:
        json.dump({
            "run_date": datetime.now().isoformat(),
            "total_samples": len(samples),
            "completed": completed,
            "failed": len(samples) - completed,
            "total_elapsed_seconds": round(time.time() - pipeline_start, 2),
            "count_matrix": matrix_path,
            "qc_json_report": json_report,
            "qc_html_report": html_report,
            "samples": all_results,
        }, f, indent=2)

    logger.info(f"Pipeline complete — {completed}/{len(samples)} samples successful")
    logger.info(f"Count matrix: {matrix_path}")
    logger.info(f"QC report: {html_report}")


if __name__ == "__main__":
    main()
