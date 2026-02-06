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

| Stage | Tool(s) | Function | Status |
|---|---|---|---|
| **0. QC** | **Fastp** | Quality filtering (Q20, 4bp sliding window). | ✅ Production |
| **1. Decontamination** | **Hostile** + **Minimap2** | Removes human/non-target DNA with aggressive short-read mapping. | ⚠️ Testing Mode* |
| **2. Classification** | **Kraken2** | Strict k-mer filtering for *Vibrio cholerae* (custom & standard DBs). | ✅ Production |
| **3. Rescue** | **MMseqs2** | **"Read Rescue"**: Recovers mutated reads via protein alignment. | ✅ Production |
| **4. Alignment** | **BWA** + **Minimap2** | Maps cleaned reads to reference genome. | ✅ Production |
| **5. Assembly** | **SPAdes** | Denovo assembly for SXT/ICE elements only (not full genome). | ⚠️ Limited Use* |
| **6. Consensus** | **Samtools** + **Pilon** | Generates consensus genome with pileup-based calling. | ✅ Production |
| **7. Polishing** | **Medaka** + **Polypolish** + **Pilon** | Platform-aware polishing: Medaka (Nanopore) or Polypolish→Pilon (Illumina). | ✅ Production |
| **8. Alignment (Final)** | **MAFFT** | Multiple sequence alignment for phylogenetic and surveillance context. | ✅ Production |
| **9. Variant Calling** | **BCFtools** | SNP and indel detection; strain differentiation with surveillance loci filtering. | ✅ Production |
| **10. Annotation** | **SnpEff** + **BLAST** | Functional annotation of variants and amr gene detection. | ✅ Production |

**Implementation Notes:**
- *Stage 1: Currently in passthrough mode for performance testing; activate for production deployment.
- *Stage 5: SPAdes is used selectively for SXT/ICE element assembly validation, not as a full denovo assembly fallback.
- *Stage 7: FMLRC2 listed in environment but not actively used; Polypolish→Pilon provides robust Illumina polishing.
- *Stage 9: Uses BCFtools mpileup/call with surveillance-aware filtering; Freebayes not currently implemented.

### 2. Intelligence Tiers (The "Brain")
*From consensus genome to actionable alert.*

| Tier | Method | Function | Deployment |
|---|---|---|---|
| **Tier 0** | **Sourmash** | **Flash Triage:** Instant k-mer drift check vs 2010/2022 baselines. | Local (CPU) |
| **Tier 1** | **HyenaDNA** | **Local AI:** Structural anomaly detection (CPU-optimized, runs locally). | Local (CPU) |
| **Tier 2** | **Evo2** (Cloud) | **Deep Forensic:** *Optional* escalation for high-risk variants. | Cloud API* |

**Cloud API Note:** Tier 2 requires Evo2 API credentials and may incur costs. Configure fallback behavior in `workflow/config/config.yaml` if cloud access is limited.

### 3. Supporting Tools

| Category | Tools | Purpose |
|---|---|---|
| **AMR Detection** | **RGI** (Resistance Gene Identifier) | Detects antibiotic resistance genes and phenotypes via sequence homology. Fallback to targeted k-mer scanning if RGI unavailable. |
| **Phylogenetics** | **IQ-TREE** + **TreeTime** | Phylogenetic tree construction (maximum likelihood) and temporal molecular clock analysis. |
| **Visualization** | **Bioconductor** (R: ggtree, treeio, ape, ggplot2) + **Matplotlib** | Tree visualization and publication-ready figures. |
| **Data Transport** | **Aiohttp** + **Httpx** | Async API calls for cloud functions and remote data retrieval. |
| **Utilities** | **Biopython** + **Pysam** + **NumPy** | Sequence parsing, SAM/BAM interface, numerical computing. |

### 4. Development Notes & Evaluated Tools

#### ⚠️ Current Development Status

**Pipeline Maturity:**
- **Production-Ready Components:** Stages 0-4, 6-10, Intelligence Tiers 0-1
- **Testing/Limited Use:** Stage 1 (Decontamination in passthrough mode), Stage 5 (Assembly for SXT elements only)
- **Cloud-Dependent:** Tier 2 (Evo2 requires API access)

**Known Limitations:**
- Hostile decontamination currently in testing mode (passthrough enabled)
- Full denovo assembly not implemented; SPAdes used only for SXT/ICE element validation
- FMLRC2 tool installed but not actively used in polishing pipeline
- Tier 2 escalation requires external API credentials and network access

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

## 📊 Pipeline Assessment

### Strengths
- **✅ Adaptive Memory Management**: Sophisticated resource tiering (BALANCED/BUNKER/EMERGENCY) enables deployment in resource-constrained settings
- **✅ Multi-Reference Strategy**: Dynamic reference selection with Haiti-specific prioritization and global fallback
- **✅ Read Rescue Innovation**: MMseqs2 protein alignment recovers mutated reads that k-mer classifiers miss
- **✅ Platform-Aware Processing**: Automatic Nanopore vs Illumina detection with appropriate polishing strategies
- **✅ Tiered Triage System**: Cost-efficient cascade (Tier 0→1→2) with early decision gates
- **✅ Transparent Documentation**: Clearly lists tools evaluated but not used, with rationale

### Current Limitations
- **⚠️ Decontamination Testing Mode**: Hostile currently in passthrough; needs activation for clinical deployment
- **⚠️ Limited Assembly Scope**: SPAdes only used for SXT elements, not full genome denovo assembly
- **⚠️ Cloud Dependency**: Tier 2 requires external API; no offline fallback for deep forensics
- **⚠️ Tool Inventory Accuracy**: Some listed tools (FMLRC2, Freebayes) not actively used in current pipeline

### Recommended for Production Deployment?
**Partial Yes, with modifications:**
- Activate Hostile decontamination before clinical use
- Verify Evo2 API configuration or implement offline Tier 2 fallback
- Consider adding full denovo assembly fallback for very low coverage samples
- Validate pipeline with representative clinical samples from target region

---

## ⚙️ Production Deployment Checklist

Before deploying this pipeline for clinical surveillance, ensure:

1. **Enable Hostile Decontamination**: The current release has decontamination in testing/passthrough mode. To activate:
   - Edit `workflow/Snakefile` line 273
   - Remove the passthrough logic and enable full decontamination
   - Verify human DNA removal with test samples

2. **Verify Evo2 API Configuration**: Tier 2 escalation requires cloud API credentials
   - Configure API keys in `workflow/config/config.yaml`
   - Note: API calls may incur costs for high-risk variant analysis
   - Consider fallback strategies if API is unavailable

3. **Database Setup**: Run `scripts/setup_databases.sh` to download:
   - Kraken2 custom Haiti database (~4GB)
   - Kraken2 standard database (~8GB)
   - MMseqs2 protein reference database (~2GB)

4. **Memory Tier Configuration**: The pipeline includes adaptive memory management
   - Default: BALANCED mode (16GB RAM)
   - BUNKER mode: <8GB RAM (reduced reference sets)
   - EMERGENCY mode: <4GB RAM (core features only)
   - Configure in `workflow/config/config.yaml`

5. **Platform Detection**: Automatic Nanopore vs Illumina detection
   - Medaka polishing for Nanopore data
   - Polypolish→Pilon pipeline for Illumina data
   - Verify correct platform detection with your sequencing setup

---

## 🙏 Acknowledgments

This pipeline builds upon and acknowledges the **[CholeraSeq pipeline](https://ceri-krisp.github.io/CholeraSeq/installation.html)**, which has contributed foundational concepts and methodologies to cholera genomic surveillance. We recognize the open-source community's efforts in advancing tools and standards for pathogen genomics.

---

## 📄 License

MIT License - See [LICENSE](LICENSE) for details.
