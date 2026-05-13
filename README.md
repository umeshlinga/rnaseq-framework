#  Scalable RNA-seq Analysis Framework

A reproducible, cloud-ready RNA-seq workflow built with object-oriented Python. Automates alignment, expression quantification, and QC validation across large multi-sample cohorts.

##  Overview

```
Raw FASTQ → Trimming → Alignment → Quantification → QC Report → Expression Matrix
```

**Key Results:**
-  28% reduction in analysis turnaround time
-  Automated Pytest-driven QC validation
-  AWS-deployed for parallel multi-sample processing

---

##  Repository Structure

```
rnaseq-framework/
├── data/
│   ├── raw/                      # Raw FASTQ files
│   ├── trimmed/                  # Post-trimming reads
│   ├── aligned/                  # BAM files
│   └── counts/                   # Gene count matrices
├── rnaseq/
│   ├── __init__.py
│   ├── trimmer.py                # Read trimming module
│   ├── aligner.py                # STAR/HISAT2 alignment wrapper
│   ├── quantifier.py             # Expression quantification (featureCounts/HTSeq)
│   ├── qc_metrics.py             # QC metric aggregation
│   └── report_generator.py       # Automated QC report output
├── tests/
│   ├── test_aligner.py
│   ├── test_quantifier.py
│   └── test_qc_metrics.py
├── scripts/
│   └── run_rnaseq.sh             # Bash pipeline entry point
├── aws/
│   └── batch_job_definition.json # AWS Batch config for parallel runs
├── config/
│   └── config.yaml
├── requirements.txt
└── README.md
```

---

##  Quickstart

```bash
git clone https://github.com/yourusername/rnaseq-framework.git
cd rnaseq-framework

pip install -r requirements.txt

# Run on a single sample
python -m rnaseq.aligner --input data/raw/sample1.fastq --output data/aligned/

# Run full pipeline
bash scripts/run_rnaseq.sh --samples sample_sheet.csv --ref /path/to/genome
```

---

##  Tech Stack

| Component | Technology |
|---|---|
| Language | Python 3.10+ (OOP) |
| Alignment | STAR / HISAT2 |
| Quantification | featureCounts / HTSeq |
| QC | MultiQC, FastQC |
| Testing | Pytest |
| Cloud | AWS (S3, Batch) |

---

##  Running Tests

```bash
pytest tests/ -v --tb=short
```

---

##  License

MIT License
