"""
rna_aligner.py
--------------
RNA-seq Alignment Module

Wraps STAR for splice-aware alignment of RNA-seq reads.
Handles genome index generation, alignment, and BAM post-processing.
"""

import os
import subprocess
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


@dataclass
class AlignmentResult:
    """Stores STAR alignment output and mapping statistics."""
    sample_id: str
    bam_path: str
    log_path: str
    total_reads: int = 0
    uniquely_mapped: int = 0
    multi_mapped: int = 0
    unmapped: int = 0
    mapping_rate: float = 0.0
    success: bool = False
    error_message: str = ""

    @property
    def summary(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "bam_path": self.bam_path,
            "total_reads": self.total_reads,
            "uniquely_mapped": self.uniquely_mapped,
            "multi_mapped": self.multi_mapped,
            "unmapped": self.unmapped,
            "mapping_rate_pct": round(self.mapping_rate * 100, 2),
            "success": self.success,
        }


class RNASeqAligner:
    """
    Aligns RNA-seq reads to a reference genome using STAR.

    STAR is preferred for RNA-seq because it performs splice-aware
    alignment, correctly handling reads spanning exon-exon junctions.

    Parameters
    ----------
    genome_dir : str
        Path to STAR genome index directory.
    output_dir : str
        Directory to write alignment output files.
    threads : int
        Number of CPU threads (default: 8).
    two_pass : bool
        Enable STAR 2-pass mode for improved novel junction detection (default: True).
    """

    def __init__(
        self,
        genome_dir: str,
        output_dir: str = "data/aligned",
        threads: int = 8,
        two_pass: bool = True,
    ):
        self.genome_dir = genome_dir
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.threads = threads
        self.two_pass = two_pass

    def _run(self, cmd: list, step: str) -> subprocess.CompletedProcess:
        logger.info(f"Running {step}")
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"{step} failed:\n{result.stderr}")
        return result

    def _parse_star_log(self, log_path: str) -> dict:
        """Parse STAR final log for mapping statistics."""
        stats = {
            "total": 0,
            "uniquely_mapped": 0,
            "multi_mapped": 0,
            "unmapped": 0,
            "mapping_rate": 0.0,
        }
        if not os.path.exists(log_path):
            return stats

        with open(log_path) as f:
            for line in f:
                line = line.strip()
                if "Number of input reads" in line:
                    stats["total"] = int(line.split("|")[-1].strip())
                elif "Uniquely mapped reads number" in line:
                    stats["uniquely_mapped"] = int(line.split("|")[-1].strip())
                elif "Uniquely mapped reads %" in line:
                    val = line.split("|")[-1].strip().replace("%", "")
                    stats["mapping_rate"] = float(val) / 100
                elif "Number of reads mapped to multiple loci" in line:
                    stats["multi_mapped"] = int(line.split("|")[-1].strip())
        stats["unmapped"] = (
            stats["total"] - stats["uniquely_mapped"] - stats["multi_mapped"]
        )
        return stats

    def build_genome_index(
        self,
        fasta: str,
        gtf: str,
        genome_dir: Optional[str] = None,
        overhang: int = 99,
    ) -> str:
        """
        Build a STAR genome index.

        Parameters
        ----------
        fasta : str
            Path to reference genome FASTA.
        gtf : str
            Path to gene annotation GTF file.
        genome_dir : str, optional
            Output directory for index (defaults to self.genome_dir).
        overhang : int
            Read length - 1. For 100bp reads, use 99.

        Returns
        -------
        str
            Path to genome index directory.
        """
        out_dir = genome_dir or self.genome_dir
        Path(out_dir).mkdir(parents=True, exist_ok=True)

        cmd = [
            "STAR",
            "--runMode", "genomeGenerate",
            "--runThreadN", str(self.threads),
            "--genomeDir", out_dir,
            "--genomeFastaFiles", fasta,
            "--sjdbGTFfile", gtf,
            "--sjdbOverhang", str(overhang),
        ]
        self._run(cmd, "STAR genome indexing")
        logger.info(f"Genome index built: {out_dir}")
        return out_dir

    def align(
        self,
        fastq_r1: str,
        sample_id: str,
        fastq_r2: Optional[str] = None,
        gtf: Optional[str] = None,
    ) -> AlignmentResult:
        """
        Align RNA-seq reads for a single sample using STAR.

        Parameters
        ----------
        fastq_r1 : str
            Path to R1 FASTQ file (or single-end FASTQ).
        sample_id : str
            Sample identifier for output naming.
        fastq_r2 : str, optional
            Path to R2 FASTQ file (paired-end).
        gtf : str, optional
            Path to GTF annotation for splice junction database.

        Returns
        -------
        AlignmentResult
        """
        out_prefix = str(self.output_dir / sample_id) + "/"
        Path(out_prefix).mkdir(parents=True, exist_ok=True)

        bam_path = out_prefix + "Aligned.sortedByCoord.out.bam"
        log_path = out_prefix + "Log.final.out"

        result = AlignmentResult(
            sample_id=sample_id,
            bam_path=bam_path,
            log_path=log_path,
        )

        read_files = [fastq_r1]
        if fastq_r2:
            read_files.append(fastq_r2)

        cmd = [
            "STAR",
            "--runThreadN", str(self.threads),
            "--genomeDir", self.genome_dir,
            "--readFilesIn", *read_files,
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--outSAMattributes", "NH", "HI", "AS", "NM", "MD",
            "--outFileNamePrefix", out_prefix,
            "--outSAMstrandField", "intronMotif",
            "--outFilterIntronMotifs", "RemoveNoncanonical",
            "--quantMode", "GeneCounts",
        ]

        # Gzipped input
        if fastq_r1.endswith(".gz"):
            cmd += ["--readFilesCommand", "zcat"]

        # GTF for novel junction detection
        if gtf:
            cmd += ["--sjdbGTFfile", gtf]

        # 2-pass mode for improved alignment
        if self.two_pass:
            cmd += ["--twopassMode", "Basic"]

        try:
            self._run(cmd, f"STAR alignment for {sample_id}")

            # Index BAM
            self._run(
                ["samtools", "index", bam_path],
                "samtools index"
            )

            # Parse alignment stats
            stats = self._parse_star_log(log_path)
            result.total_reads = stats["total"]
            result.uniquely_mapped = stats["uniquely_mapped"]
            result.multi_mapped = stats["multi_mapped"]
            result.unmapped = stats["unmapped"]
            result.mapping_rate = stats["mapping_rate"]
            result.success = True

            logger.info(
                f"Alignment complete — {result.uniquely_mapped} uniquely mapped reads "
                f"({result.mapping_rate*100:.1f}%)"
            )

        except Exception as e:
            result.error_message = str(e)
            logger.error(f"Alignment failed for {sample_id}: {e}")

        return result
