# Vibrion Sentinel

**Public Health Decision Support System for Cholera Surveillance**
**Version:** 2.0

Vibrion Sentinel is an advanced genomic surveillance pipeline designed for the identification, characterization, and molecular epidemiology of *Vibrio cholerae*. Originally developed to support response efforts during the 2022 Haiti Resurgence, Version 2.0 introduces reference-agnostic capabilities and "Salvage Mode" to maximize data recovery from complex or low-quality field samples.

## Key Features (Version 2.0)

### 1. Salvage Mode for Low-Quality Data
Standard pipelines often discard short reads (<50bp) or unclassified sequences. Vibrion Sentinel integrates **Evo2** and **MMseqs2** to rescue these reads, enabling successful consensus generation from degraded samples or those with high host contamination. It includes robust logic for paired-end file handling and broken mate-pair recovery.

### 2. Public Health Typing
The pipeline automatically extracts and interprets critical markers for public health decision-making:
*   **Serotype Prediction:** Analyzes *wbeT* frameshifts to distinguish Ogawa (wild-type) from Inaba (mutant) serotypes.
*   **Toxin Genotyping:** Identifies *ctxB* alleles (e.g., ctxB7 vs. ctxB1) to assess virulence potential.
*   **Antimicrobial Resistance (AMR):** Scans for resistance determinants in the SXT element and plasmids to guide antibiotic stewardship (e.g., Doxycycline susceptibility).

### 3. Heterogeneity Detection
To support complex outbreak analysis, the pipeline identifies sub-clonal populations (minor variants with allele frequency 10-90%). This functionality is critical for detecting:
*   Mixed infections (co-infection with multiple strains).
*   Rapid in-host evolution (e.g., active serotype switching).
*   Sample contamination.

### 4. Forensic Validation
Ensures data integrity through:
*   **Housekeeping Checksum:** Validates assembly integrity against 7PET markers (*recA*, *gyrB*, *dnaE*).
*   **Phylogenetic Placement:** Automatically places samples on a global phylogenetic tree to distinguish local resurgence from foreign importation.

## Installation

### Prerequisites
*   **Operating System:** Linux or macOS
*   **Package Manager:** Conda or Mamba
*   **Hardware:** Minimum 16GB RAM recommended (8GB minimum).

### Steps
1.  **Clone the repository:**
    ```bash
    git clone https://github.com/intelogroup/vibrion-sentinel.git
    cd vibrion-sentinel
    ```

2.  **Create the environment:**
    ```bash
    conda env create -f environment.yml
    conda activate vibrion
    ```

3.  **Prepare References:**
    Ensure reference databases (Kraken2, MMseqs2) and genome assemblies (2010EL-1786) are available in the `data/` directory.

## Configuration

### Environment Variables
For the Evo2 AI rescue feature to function, you must provide an NVIDIA API key. For security, **do not commit this key**. Set it as an environment variable before running the pipeline:

```bash
export NVIDIA_API_KEY="nvapi-..."
```

### Pipeline Settings
Edit `workflow/config/config.yaml` to adjust run parameters:

```yaml
# Directory Paths
samples_dir: "data/raw_reads"
output_dir: "data/pipeline_output"

# Analysis Thresholds
mapping_rescue_enabled: True
mapping_rescue_confidence_threshold: 0.6
```

## Usage

The pipeline is built on Snakemake.

**Standard Run (4 cores):**
```bash
snakemake --use-conda --cores 4
```

**Run for a specific sample:**
```bash
snakemake --cores 4 data/pipeline_output/SRR24010030/08_comprehensive_report/surveillance_report.md
```

## Outputs

Results are organized by sample ID in `data/pipeline_output/{SAMPLE_ID}/`:

*   **`08_comprehensive_report/`**: Contains `surveillance_report.md`, the primary actionable report.
*   **`09_consensus/`**: The polished consensus genome (*.fasta) and assembly metrics.
*   **`05_variants/`**: VCF files and SNP reports detailing variants found against the reference.
*   **`06_amr/`**: Antibiotic resistance profiles and virulence factors.
*   **`10_phylogeny/`**: Phylogenetic tree image showing the sample's placement in the global context.

## License

MIT License. See `LICENSE` file for details.