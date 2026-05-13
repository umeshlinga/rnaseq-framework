"""
quantifier.py
-------------
Gene Expression Quantification Module

Wraps featureCounts for read counting across genomic features.
Produces a gene-level count matrix suitable for differential expression analysis.
"""

import subprocess
import logging
import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Dict, Optional

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


@dataclass
class QuantificationResult:
    """Stores quantification output and statistics."""
    sample_id: str
    counts_path: str
    total_fragments: int = 0
    assigned_fragments: int = 0
    unassigned_multimapping: int = 0
    unassigned_no_features: int = 0
    assignment_rate: float = 0.0
    success: bool = False
    error_message: str = ""

    @property
    def summary(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "counts_path": self.counts_path,
            "total_fragments": self.total_fragments,
            "assigned_fragments": self.assigned_fragments,
            "assignment_rate_pct": round(self.assignment_rate * 100, 2),
            "unassigned_multimapping": self.unassigned_multimapping,
            "unassigned_no_features": self.unassigned_no_features,
            "success": self.success,
        }


class GeneQuantifier:
    """
    Quantifies gene expression from BAM files using featureCounts.

    featureCounts assigns aligned reads to genomic features (genes/exons)
    defined in a GTF annotation file.

    Parameters
    ----------
    gtf : str
        Path to gene annotation GTF file.
    output_dir : str
        Directory to write count files.
    feature_type : str
        Feature type to quantify (default: 'exon').
    attribute_type : str
        GTF attribute used as feature ID (default: 'gene_id').
    paired_end : bool
        Whether reads are paired-end (default: True).
    strand_specific : int
        Strand specificity: 0=unstranded, 1=stranded, 2=reverse-stranded.
    threads : int
        Number of threads (default: 4).
    min_mapping_quality : int
        Minimum mapping quality for read assignment (default: 10).
    """

    def __init__(
        self,
        gtf: str,
        output_dir: str = "data/counts",
        feature_type: str = "exon",
        attribute_type: str = "gene_id",
        paired_end: bool = True,
        strand_specific: int = 2,
        threads: int = 4,
        min_mapping_quality: int = 10,
    ):
        self.gtf = gtf
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.feature_type = feature_type
        self.attribute_type = attribute_type
        self.paired_end = paired_end
        self.strand_specific = strand_specific
        self.threads = threads
        self.min_mapping_quality = min_mapping_quality

    def _parse_summary(self, summary_path: str) -> dict:
        """Parse featureCounts summary file for assignment statistics."""
        stats = {
            "assigned": 0,
            "unassigned_multimapping": 0,
            "unassigned_nofeatures": 0,
            "total": 0,
        }
        if not Path(summary_path).exists():
            return stats

        with open(summary_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                status = row.get("Status", "")
                count = int(list(row.values())[-1]) if len(row) > 1 else 0
                if status == "Assigned":
                    stats["assigned"] = count
                elif status == "Unassigned_MultiMapping":
                    stats["unassigned_multimapping"] = count
                elif status == "Unassigned_NoFeatures":
                    stats["unassigned_nofeatures"] = count
                stats["total"] += count
        return stats

    def quantify(self, bam_path: str, sample_id: str) -> QuantificationResult:
        """
        Count reads per gene for a single sample.

        Parameters
        ----------
        bam_path : str
            Path to sorted, indexed BAM file.
        sample_id : str
            Sample identifier for output naming.

        Returns
        -------
        QuantificationResult
        """
        counts_path = str(self.output_dir / f"{sample_id}.counts.txt")
        summary_path = counts_path + ".summary"

        result = QuantificationResult(
            sample_id=sample_id,
            counts_path=counts_path,
        )

        cmd = [
            "featureCounts",
            "-T", str(self.threads),
            "-t", self.feature_type,
            "-g", self.attribute_type,
            "-s", str(self.strand_specific),
            "-Q", str(self.min_mapping_quality),
            "-a", self.gtf,
            "-o", counts_path,
        ]

        if self.paired_end:
            cmd += ["-p", "--countReadPairs"]

        cmd.append(bam_path)

        logger.info(f"Quantifying gene expression for {sample_id}")
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if proc.returncode != 0:
            result.error_message = proc.stderr
            logger.error(f"featureCounts failed for {sample_id}: {proc.stderr}")
            return result

        stats = self._parse_summary(summary_path)
        result.total_fragments = stats["total"]
        result.assigned_fragments = stats["assigned"]
        result.unassigned_multimapping = stats["unassigned_multimapping"]
        result.unassigned_no_features = stats["unassigned_nofeatures"]
        if result.total_fragments > 0:
            result.assignment_rate = result.assigned_fragments / result.total_fragments
        result.success = True

        logger.info(
            f"Quantification complete — {result.assigned_fragments} assigned fragments "
            f"({result.assignment_rate*100:.1f}%)"
        )
        return result

    def build_count_matrix(
        self,
        count_files: List[str],
        sample_ids: List[str],
        output_path: Optional[str] = None,
    ) -> Dict[str, Dict[str, int]]:
        """
        Merge per-sample count files into a unified count matrix.

        Parameters
        ----------
        count_files : list of str
            Paths to per-sample featureCounts output files.
        sample_ids : list of str
            Sample identifiers corresponding to each count file.
        output_path : str, optional
            If provided, write the matrix to this TSV file.

        Returns
        -------
        dict
            Nested dict: {gene_id: {sample_id: count}}.
        """
        matrix: Dict[str, Dict[str, int]] = {}

        for counts_file, sample_id in zip(count_files, sample_ids):
            with open(counts_file) as f:
                for line in f:
                    if line.startswith("#") or line.startswith("Geneid"):
                        continue
                    parts = line.strip().split("\t")
                    gene_id = parts[0]
                    count = int(parts[-1])
                    if gene_id not in matrix:
                        matrix[gene_id] = {}
                    matrix[gene_id][sample_id] = count

        if output_path:
            with open(output_path, "w") as out:
                header = "\t".join(["gene_id"] + sample_ids)
                out.write(header + "\n")
                for gene_id, counts in sorted(matrix.items()):
                    row = [gene_id] + [str(counts.get(s, 0)) for s in sample_ids]
                    out.write("\t".join(row) + "\n")
            logger.info(f"Count matrix written: {output_path}")

        logger.info(
            f"Count matrix built — {len(matrix)} genes x {len(sample_ids)} samples"
        )
        return matrix

    def compute_rpkm(
        self,
        count_matrix: Dict[str, Dict[str, int]],
        gene_lengths: Dict[str, int],
    ) -> Dict[str, Dict[str, float]]:
        """
        Compute RPKM (Reads Per Kilobase per Million mapped reads).

        Parameters
        ----------
        count_matrix : dict
            Gene x sample count matrix.
        gene_lengths : dict
            Gene lengths in base pairs {gene_id: length}.

        Returns
        -------
        dict
            RPKM values: {gene_id: {sample_id: rpkm}}.
        """
        rpkm: Dict[str, Dict[str, float]] = {}

        # Compute per-sample total mapped reads (millions)
        sample_totals: Dict[str, float] = {}
        all_samples = set()
        for gene_counts in count_matrix.values():
            all_samples.update(gene_counts.keys())
        for sample in all_samples:
            sample_totals[sample] = sum(
                counts.get(sample, 0) for counts in count_matrix.values()
            ) / 1e6

        for gene_id, gene_counts in count_matrix.items():
            length_kb = gene_lengths.get(gene_id, 1000) / 1000
            rpkm[gene_id] = {}
            for sample_id, count in gene_counts.items():
                total_m = sample_totals.get(sample_id, 1.0)
                rpkm[gene_id][sample_id] = round(count / (length_kb * total_m), 4)

        return rpkm
