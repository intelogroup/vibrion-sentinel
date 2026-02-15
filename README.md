# Vibrion Sentinel 🛡️🇭🇹

**Public Health Decision Support System for Cholera Surveillance**
**Version:** 2.0 (Agentic Sentinel)

Vibrion Sentinel is an advanced genomic surveillance pipeline designed for the rapid identification and characterization of *Vibrio cholerae* in outbreak settings. Originally developed for the 2022 Haiti Resurgence, it has evolved into a reference-agnostic system capable of handling complex, low-quality, and mixed-infection samples.

## 🚀 Key Features (v2.0)

### 1. 🧬 Salvage Mode (New!)
**AI-Powered Rescue for Low-Quality Data**
- **Evo2 & MMseqs2 Integration:** Automatically recovers unclassified or short reads rejected by standard classifiers.
- **Robust Paired-End Handling:** Intelligent pairing logic salvages data even when mate pairs are corrupted or missing.
- **Zero-Loss Consensus:** Generates alignment-based consensus even when *de novo* assembly fails due to low coverage.

### 2. 🏥 Public Health Typing
**Clinical Decision Support**
- **Serotype Prediction:** Detects `wbeT` frameshifts to distinguish Ogawa (wild-type) from Inaba (mutant) serotypes, critical for vaccine deployment.
- **Toxin Genotyping:** Identifies `ctxB` alleles (e.g., `ctxB7` Haiti/Classical vs `ctxB1` El Tor) to assess virulence potential.
- **AMR Profiling:** Scans for resistance genes (SXT element, plasmids) to guide antibiotic treatment (Doxycycline/Azithromycin/Ciprofloxacin).

### 3. 🔍 Heterogeneity Detection
**Mixed Infection Alert System**
- **Minor Variant Calling:** Identifies sub-clonal populations (AF 10-90%) indicating mixed infections or rapid in-host evolution.
- **Strain Triage:** Distinguishes between outbreak strains and environmental non-O1/O139 lineages.

### 4. 🛡️ Forensic Validation
- **Housekeeping Checksum:** Validates assembly integrity against 7PET markers (`recA`, `gyrB`, `dnaE`).
- **Phylogenetic Placement:** Places samples on a global tree to identify origin (e.g., local resurgence vs. foreign import).

---

## 📦 Installation

### Prerequisites
- **Conda/Mamba:** Recommended for environment management.
- **NVIDIA API Key:** Required for Evo2 AI rescue (optional but recommended).

### 1. Clone the Repository
```bash
git clone https://github.com/intelogroup/vibrion-sentinel.git
cd vibrion-sentinel
```

### 2. Create Environment
```bash
conda env create -f environment.yml
conda activate vibrion
```

### 3. Configure References
Ensure reference databases are available. You may need to build or download them:
- **Kraken2 Database:** Standard or Custom Haiti build.
- **MMseqs2 Database:** SwissProt or similar.
- **Reference Genome:** *V. cholerae* 2010EL-1786 (CP003069.1/CP003070.1).

---

## ⚙️ Configuration & Security

**⚠️ SECURITY WARNING: NEVER commit API keys to version control.**

The pipeline requires an NVIDIA API key for the Evo2 rescue feature. Set this as an environment variable before running the pipeline:

```bash
export NVIDIA_API_KEY="nvapi-..."
```

Alternatively, you can pass it via configuration file, but ensure that file is **excluded** from git (add to `.gitignore`).

### Config File (`workflow/config/config.yaml`)
Adjust paths and thresholds in the configuration file:

```yaml
# Input/Output
samples_dir: "data/raw_reads"
output_dir: "data/pipeline_output"

# Thresholds
mapping_rescue_enabled: True
mapping_rescue_confidence_threshold: 0.6
```

---

## 🏃‍♂️ Usage

**Run the pipeline with Snakemake:**

```bash
# Run locally with 4 cores
snakemake --use-conda --cores 4

# Run specific sample
snakemake --cores 4 data/pipeline_output/SRR24010030/08_comprehensive_report/surveillance_report.md
```

### Output Artifacts
Results are stored in `data/pipeline_output/{SAMPLE_ID}/`:

- **`08_comprehensive_report/surveillance_report.md`**: The primary actionable report for public health officials.
- **`09_consensus/`**: Polished consensus genomes (`*_polished.fasta`).
- **`05_variants/`**: VCF files and SNP reports.
- **`06_amr/`**: Antibiotic resistance profiles.
- **`10_phylogeny/`**: Phylogenetic tree placement (`tree.png`).

---

## 🧪 Testing & Validation

This repository includes a stress test framework to validate pipeline integrity.

**Run Validation Harness:**
```bash
python3 validation/validation_harness.py
```

See `validation/STRESS_TEST_PROTOCOL.md` for detailed objectives.

---

## 📜 License

MIT License. See `LICENSE` for details.

---

**Repository:** https://github.com/intelogroup/vibrion-sentinel
**Status:** ✅ Field-Ready (v2.0)
**Maintainer:** Intelo Group
