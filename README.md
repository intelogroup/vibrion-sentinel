# 🛡️ Vibrion Sentinel v2.0 (Public Release)

**Clinical-Grade Genomic Surveillance for Cholera in Resource-Constrained Settings**

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

> *"Read Rescue" + "Local-First Philosophy" + "Forensic Gold Standard"*

---

## 🎯 Mission

Vibrion Sentinel provides **real-time cholera strain surveillance**, designed for:
- 🏥 **Clinical decision support** (stop-light reporting)

This repository contains the **Production Pipeline** code, stripped of internal R&D scripts, ready for deployment.

---

## 📦 Installation

### Prerequisites
- **Miniforge3** (Mamba/Conda)
- **Docker** (Optional, for containerized run)
- 16GB RAM minimum

### 1. Clone & Setup
```bash
git clone https://github.com/your-org/vibrion-sentinel.git
cd vibrion-sentinel

# Create Conda Environment
mamba env create -f environment.yml
conda activate vibrion-sentinel
```

### 2. Download Databases
The pipeline requires specific reference databases (~10GB).
```bash
# Download and build Kraken2 & MMseqs2 databases
bash scripts/setup_databases.sh
```

---

## 🚀 Usage

### Place your samples
Put your `.fastq.gz` sequencing files in `data/raw_reads/`.

### Run the Pipeline
```bash
# Analyze a specific sample
snakemake --cores 4 --use-conda \
  --config samples='["YOUR_SAMPLE_ID"]' \
  data/pipeline_output/YOUR_SAMPLE_ID/08_comprehensive_report/surveillance_report.md
```

### Output
Results will be in `data/pipeline_output/YOUR_SAMPLE_ID/`.
- **Report:** `08_comprehensive_report/surveillance_report.md` (Forensic Summary)
- **Tree:** `10_phylogeny/tree.png`
- **Triage:** `07_triage/triage_decision.json`

---

## 🛡️ System Architecture

### 1. Data Pipeline (The "Body")
*From raw dirty water to a clean forensic genome.*

| Stage | Tool | Function |
|---|---|---|
| **1. Decon** | Hostile | Removes human/non-target DNA. |
| **2. Classify** | Kraken2 | Strict k-mer filtering for *Vibrio cholerae*. |
| **3. Rescue** | **MMseqs2** | **"Read Rescue"**: Recovers mutated reads via protein alignment. |
| **4. Assembly** | BWA/Pilon | Builds the consensus genome for analysis. |

### 2. Intelligence Tiers (The "Brain")
*From consensus genome to actionable alert.*

| Tier | Method | Function |
|---|---|---|
| **Tier 0** | Sourmash | **Flash Triage:** Instant k-mer drift check vs 2010/2022 baselines. |
| **Tier 1** | HyenaDNA | **Local AI:** Structural anomaly detection (CPU-optimized). |
| **Tier 2** | Evo2 (Cloud) | **Deep Forensic:** *Optional* escalation for high-risk variants. |

---

## 📄 License

MIT License - See [LICENSE](LICENSE) for details.
