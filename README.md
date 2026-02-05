<div align="center">
  
  # <img src="assets/logo.png" alt="Vibrion Sentinel Logo" width="32" style="margin-right: 10px; vertical-align: middle;" /> Vibrion Sentinel v2.0 (Public Release)
  
  **Clinical-Grade Genomic Surveillance for Cholera in Resource-Constrained Settings**
  
  *Built upon the [CholeraSeq pipeline](https://ceri-krisp.github.io/CholeraSeq/installation.html) — advancing open-source cholera genomics.*
  
  [![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
  
  > *"Read Rescue" + "Local-First Philosophy" + "Forensic Gold Standard"*
</div>

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

| Stage | Tool(s) | Function |
|---|---|---|
| **0. QC** | **Fastp** | Quality filtering (Q20, 4bp sliding window). |
| **1. Decontamination** | **Hostile** + **Minimap2** | Removes human/non-target DNA with aggressive short-read mapping. |
| **2. Classification** | **Kraken2** | Strict k-mer filtering for *Vibrio cholerae* (custom & standard DBs). |
| **3. Rescue** | **MMseqs2** | **"Read Rescue"**: Recovers mutated reads via protein alignment. |
| **4. Alignment** | **BWA** + **Minimap2** | Maps cleaned reads to reference genome. |
| **5. Assembly** | **SPAdes** | Denovo assembly backup if coverage too low. |
| **6. Consensus** | **Samtools** + **Pilon** | Generates consensus genome with pileup-based calling. |
| **7. Polishing** | **Medaka** + **Polypolish** + **FMLRC2** | Iterative error correction; Nanopore-aware (medaka) or Illumina-aware (pilon) polishing. |
| **8. Alignment (Final)** | **MAFFT** | Multiple sequence alignment for phylogenetic and surveillance context. |
| **9. Variant Calling** | **BCFtools** + **Freebayes** | SNP and indel detection; strain differentiation. |
| **10. Annotation** | **SnpEff** + **BLAST** | Functional annotation of variants and amr gene detection. |

### 2. Intelligence Tiers (The "Brain")
*From consensus genome to actionable alert.*

| Tier | Method | Function |
|---|---|---|
| **Tier 0** | **Sourmash** | **Flash Triage:** Instant k-mer drift check vs 2010/2022 baselines. |
| **Tier 1** | **HyenaDNA** | **Local AI:** Structural anomaly detection (CPU-optimized, runs locally). |
| **Tier 2** | **Evo2** (Cloud) | **Deep Forensic:** *Optional* escalation for high-risk variants. |

### 3. Supporting Tools

| Category | Tools | Purpose |
|---|---|---|
| **AMR Detection** | **RGI** (Resistance Gene Identifier) | Detects antibiotic resistance genes and phenotypes via sequence homology. Fallback to targeted k-mer scanning if RGI unavailable. |
| **Phylogenetics** | **IQ-TREE** + **TreeTime** | Phylogenetic tree construction (maximum likelihood) and temporal molecular clock analysis. |
| **Visualization** | **Bioconductor** (R: ggtree, treeio, ape, ggplot2) + **Matplotlib** | Tree visualization and publication-ready figures. |
| **Data Transport** | **Aiohttp** + **Httpx** | Async API calls for cloud functions and remote data retrieval. |
| **Utilities** | **Biopython** + **Pysam** + **NumPy** | Sequence parsing, SAM/BAM interface, numerical computing. |

### 4. Tools Evaluated But Not Used

#### ❌ **Caduceus** (Nucleotide Transformer)
- **Issue:** Missing dependency (`mamba_ssm` package not available in stable environments)
- **Intended role:** Tier 2 locus structural verification
- **Impact:** Tier 2 layer non-functional; all samples auto-escalate to Evo2
- **Decision:** Focus on proven Evo2 API rather than debug underdeveloped models

#### ❌ **Centrifuge** (Metagenomic Classifier)
- **Issue:** Over-engineered for this use case; requires separate database build/maintenance
- **Trade-offs:** 
  - Slower than Kraken2 (more granular but computationally expensive)
  - Requires redundant taxonomy building
  - Minimal sensitivity gain for cholera surveillance (specialized lineages)
- **Decision:** Kraken2 + MMseqs2 "Read Rescue" provides sufficient sensitivity with simpler architecture

---

## 🙏 Acknowledgments

This pipeline builds upon and acknowledges the **[CholeraSeq pipeline](https://ceri-krisp.github.io/CholeraSeq/installation.html)**, which has contributed foundational concepts and methodologies to cholera genomic surveillance. We recognize the open-source community's efforts in advancing tools and standards for pathogen genomics.

---

## 📄 License

MIT License - See [LICENSE](LICENSE) for details.
