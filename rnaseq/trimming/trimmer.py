"""
trimmer.py
----------
Read Trimming Module

Wraps Trimmomatic for adapter trimming and quality-based filtering.
Supports both single-end and paired-end FASTQ files.
"""

import subprocess
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


@dataclass
class TrimmingResult:
    """Stores trimming output and statistics."""
    sample_id: str
    trimmed_r1: str
    trimmed_r2: Optional[str] = None
    input_reads: int = 0
    surviving_reads: int = 0
    dropped_reads: int = 0
    success: bool = False
    error_message: str = ""

    @property
    def survival_rate(self) -> float:
        if self.input_reads == 0:
            return 0.0
        return round((self.surviving_reads / self.input_reads) * 100, 2)

    @property
    def summary(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "trimmed_r1": self.trimmed_r1,
            "trimmed_r2": self.trimmed_r2,
            "input_reads": self.input_reads,
            "surviving_reads": self.surviving_reads,
            "survival_rate_pct": self.survival_rate,
            "dropped_reads": self.dropped_reads,
            "success": self.success,
        }


class ReadTrimmer:
    """
    Trims adapter sequences and low-quality bases using Trimmomatic.

    Handles both single-end and paired-end RNA-seq FASTQ files.
    Supports gzipped input and produces gzipped output.

    Parameters
    ----------
    adapters : str
        Path to adapter FASTA file (e.g. TruSeq3-PE.fa).
    output_dir : str
        Directory to write trimmed FASTQ files.
    leading : int
        Remove leading bases below this quality (default: 3).
    trailing : int
        Remove trailing bases below this quality (default: 3).
    sliding_window : str
        Sliding window trimming: 'window_size:required_quality' (default: '4:15').
    min_length : int
        Minimum read length after trimming (default: 36).
    threads : int
        Number of threads (default: 4).
    """

    def __init__(
        self,
        adapters: str,
        output_dir: str = "data/trimmed",
        leading: int = 3,
        trailing: int = 3,
        sliding_window: str = "4:15",
        min_length: int = 36,
        threads: int = 4,
    ):
        self.adapters = adapters
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.leading = leading
        self.trailing = trailing
        self.sliding_window = sliding_window
        self.min_length = min_length
        self.threads = threads

    def _parse_trimmomatic_log(self, stderr: str) -> dict:
        """Parse Trimmomatic stderr for read survival statistics."""
        stats = {"input": 0, "surviving": 0, "dropped": 0}
        for line in stderr.splitlines():
            if "Input Read Pairs:" in line:
                parts = line.split()
                try:
                    stats["input"] = int(parts[3])
                    stats["surviving"] = int(parts[6])
                    stats["dropped"] = int(parts[17])
                except (IndexError, ValueError):
                    pass
            elif "Input Reads:" in line:
                parts = line.split()
                try:
                    stats["input"] = int(parts[2])
                    stats["surviving"] = int(parts[4])
                    stats["dropped"] = int(parts[8])
                except (IndexError, ValueError):
                    pass
        return stats

    def trim_paired(
        self,
        fastq_r1: str,
        fastq_r2: str,
        sample_id: str,
    ) -> TrimmingResult:
        """
        Trim paired-end reads.

        Parameters
        ----------
        fastq_r1 : str
            Path to R1 FASTQ file.
        fastq_r2 : str
            Path to R2 FASTQ file.
        sample_id : str
            Sample identifier for output naming.

        Returns
        -------
        TrimmingResult
        """
        out_r1 = str(self.output_dir / f"{sample_id}_R1_trimmed.fastq.gz")
        out_r2 = str(self.output_dir / f"{sample_id}_R2_trimmed.fastq.gz")
        out_r1_unpaired = str(self.output_dir / f"{sample_id}_R1_unpaired.fastq.gz")
        out_r2_unpaired = str(self.output_dir / f"{sample_id}_R2_unpaired.fastq.gz")

        result = TrimmingResult(
            sample_id=sample_id,
            trimmed_r1=out_r1,
            trimmed_r2=out_r2,
        )

        cmd = [
            "trimmomatic", "PE",
            "-threads", str(self.threads),
            "-phred33",
            fastq_r1, fastq_r2,
            out_r1, out_r1_unpaired,
            out_r2, out_r2_unpaired,
            f"ILLUMINACLIP:{self.adapters}:2:30:10:2:keepBothReads",
            f"LEADING:{self.leading}",
            f"TRAILING:{self.trailing}",
            f"SLIDINGWINDOW:{self.sliding_window}",
            f"MINLEN:{self.min_length}",
        ]

        logger.info(f"Trimming paired-end reads for {sample_id}")
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if proc.returncode != 0:
            result.error_message = proc.stderr
            logger.error(f"Trimmomatic failed for {sample_id}: {proc.stderr}")
            return result

        stats = self._parse_trimmomatic_log(proc.stderr)
        result.input_reads = stats["input"]
        result.surviving_reads = stats["surviving"]
        result.dropped_reads = stats["dropped"]
        result.success = True

        logger.info(
            f"Trimming complete — {result.surviving_reads}/{result.input_reads} "
            f"read pairs surviving ({result.survival_rate}%)"
        )
        return result

    def trim_single(self, fastq: str, sample_id: str) -> TrimmingResult:
        """
        Trim single-end reads.

        Parameters
        ----------
        fastq : str
            Path to input FASTQ file.
        sample_id : str
            Sample identifier.

        Returns
        -------
        TrimmingResult
        """
        out = str(self.output_dir / f"{sample_id}_trimmed.fastq.gz")
        result = TrimmingResult(sample_id=sample_id, trimmed_r1=out)

        cmd = [
            "trimmomatic", "SE",
            "-threads", str(self.threads),
            "-phred33",
            fastq, out,
            f"ILLUMINACLIP:{self.adapters}:2:30:10",
            f"LEADING:{self.leading}",
            f"TRAILING:{self.trailing}",
            f"SLIDINGWINDOW:{self.sliding_window}",
            f"MINLEN:{self.min_length}",
        ]

        logger.info(f"Trimming single-end reads for {sample_id}")
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if proc.returncode != 0:
            result.error_message = proc.stderr
            logger.error(f"Trimmomatic failed for {sample_id}: {proc.stderr}")
            return result

        stats = self._parse_trimmomatic_log(proc.stderr)
        result.input_reads = stats["input"]
        result.surviving_reads = stats["surviving"]
        result.dropped_reads = stats["dropped"]
        result.success = True

        logger.info(
            f"Trimming complete — {result.surviving_reads}/{result.input_reads} "
            f"reads surviving ({result.survival_rate}%)"
        )
        return result
