# Vibrion Sentinel

**Genomic Surveillance Pipeline for *Vibrio cholerae* — v2.0**

Vibrion Sentinel is a field-deployable Snakemake pipeline for cholera genomic surveillance. It accepts raw Illumina or Nanopore FASTQ reads and produces a full public-health report: serotype, toxin genotype, AMR profile, CTXφ integration status, phylogenetic placement, and a **Vibrio Resurgence Score (VRS)**.

Originally developed to support the 2022 Haiti resurgence response. v2.0 adds reference-agnostic triage, AI-assisted anomaly detection (HyenaDNA + Evo2), and a tiered memory management system for deployment on laptops down to 6GB RAM.

---

## Quick Start

```bash
# 1. Clone
git clone https://github.com/intelogroup/vibrion-sentinel.git
cd vibrion-sentinel

# 2. Create environment
conda env create -f environment.yml
conda activate vibrion

# 3. Download databases (~10GB total)
bash scripts/setup_databases.sh

# 4. Place reads
cp /path/to/your/sample.fastq.gz data/raw_reads/MY_SAMPLE.fastq.gz

# 5. Run
export NVIDIA_API_KEY="nvapi-..."   # optional — enables cloud Evo2 escalation
bash scripts/run_pipeline.sh --config workflow/test_config.yaml --cores 8
```

Results → `data/pipeline_output/MY_SAMPLE/08_comprehensive_report/surveillance_report.md`

---

## System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| RAM | 8 GB | 16 GB |
| CPU cores | 2 | 8 |
| Disk | 15 GB free | 30 GB free |
| OS | Linux, macOS | Linux, macOS |
| Conda/Mamba | required | mamba (faster) |

The pipeline auto-detects available RAM and selects a memory tier (FULL / BALANCED / BUNKER / EMERGENCY). See `workflow/test_config.yaml` for tier definitions.

---

## Installation

### 1. Conda environment

```bash
conda env create -f environment.yml
conda activate vibrion
```

Tools installed: `snakemake`, `fastp`, `bwa`, `samtools`, `kraken2`, `bcftools`, `spades`, `pilon`, `sourmash`, `mmseqs2`, `snpEff`, `mafft`, `FastTree`, `hostile`.

### 2. HyenaDNA (AI triage — Tier 1)

HyenaDNA runs locally for genomic anomaly scoring without cloud API calls. See **[docs/HYENADNA_SETUP.md](docs/HYENADNA_SETUP.md)** for the full guide.

**TL;DR:** included in `environment.yml`. Verify with:
```bash
python3 -c "from transformers import AutoTokenizer; print('HyenaDNA deps OK')"
```

To disable and fall back to k-mer-only triage:
```yaml
triage:
  hyena_use_real_model: false
```

### 3. Databases

```bash
bash scripts/setup_databases.sh
```

| Database | Size | Purpose |
|----------|------|---------|
| Kraken2 standard 8GB | 8.0 GB | Taxonomic classification |
| Kraken2 serogroup DB | ~20 MB | O1/O139 serogroup probing (built locally) |
| MMseqs2 SwissProt | ~1.0 GB | Unclassified read rescue |
| Reference genomes | ~100 MB | 2010EL-1786, Haiti 2022, global panel |
| Global references | ~2.1 GB | DRC-2024, Yemen, India Wave3, Malawi, S. Africa, Bangladesh |
| Core alignment | ~8 MB | Phylogenetic backbone |

Skip large downloads for quick testing:
```bash
bash scripts/setup_databases.sh --skip-kraken --skip-mmseqs
```

---

## NVIDIA API Key (optional)

Evo2 cloud escalation is **optional**. The pipeline runs fully locally — Evo2 is only called when Tier 0 (sourmash k-mer) + Tier 1 (HyenaDNA) triage flags a potential anomaly. Most samples never reach it.

Get a free key at **https://build.nvidia.com/arc/evo2** then:
```bash
export NVIDIA_API_KEY="nvapi-..."
```

> **Security:** Never commit your key. All configs in this repo have `nvidia_api_key: ""`.  
> The runner script injects it at runtime via `--config nvidia_api_key=$NVIDIA_API_KEY`.

---

## Usage

```bash
# Validate first (dry-run)
bash scripts/run_pipeline.sh --config workflow/test_config.yaml --dry-run

# Full run
bash scripts/run_pipeline.sh --config workflow/test_config.yaml --cores 8
```

### Writing a config

```yaml
samples_dir: "data/raw_reads"
output_dir: "data/pipeline_output"
reference_dir: "data/references"
global_references_dir: "data/global_references"
kraken_db: "data/kraken2_standard_8gb"
serogroup_db: "data/kraken2_serogroup"
threads: 8
memory_mb: 16000
pipeline_mode: "LABORATORY_FULL"
nvidia_api_key: ""   # do not hardcode — injected by run_pipeline.sh
samples:
  - MY_SAMPLE   # matches data/raw_reads/MY_SAMPLE.fastq.gz
                # or MY_SAMPLE_1.fastq.gz + MY_SAMPLE_2.fastq.gz (paired-end)
```

### Input formats

| Format | File convention |
|--------|----------------|
| Single-end | `data/raw_reads/SAMPLE.fastq.gz` |
| Paired-end | `data/raw_reads/SAMPLE_1.fastq.gz` + `SAMPLE_2.fastq.gz` |
| Nanopore | single-end `.fastq.gz` — platform auto-detected |

---

## Pipeline Overview

```
Raw FASTQ
  ├─ fastp QC
  ├─ hostile (human decontamination)
  ├─ Kraken2 classification
  ├─ Vibrio read extraction + MMseqs2 rescue of unclassified reads
  ├─ Serogroup detection (O1 / O139 / NOVC)
  ├─ BWA alignment → auto-selected regional reference
  ├─ Variant calling (bcftools) + SnpEff annotation
  ├─ CTXφ phage integration detection (dual dif-site)
  ├─ SXT element assembly
  ├─ AMR profiling (targeted + RGI)
  ├─ Phenotypic prediction (biofilm, motility, rugose state)
  ├─ Consensus genome + Pilon polishing
  ├─ Phylogenetic placement (MAFFT + FastTree)
  ├─ Triage: Tier 0 sourmash → Tier 1 HyenaDNA → Tier 2 Evo2 (cloud, if needed)
  ├─ VRS (Vibrio Resurgence Score) calculation
  └─ Comprehensive surveillance report (Markdown + JSON)
```

---

## Outputs

All outputs are in `data/pipeline_output/<SAMPLE_ID>/`:

| Directory | Key files |
|-----------|-----------|
| `08_comprehensive_report/` | `surveillance_report.md`, `vrs_score.json` |
| `09_consensus/` | `*_polished.fasta`, `ctx_integration.json`, `platform_detection.json` |
| `05_variants/` | `*.vcf.gz`, `snp_report.json`, `haplotypes.json`, `public_health_typing.json` |
| `06_amr/` | `amr_report.json`, `rgi_report.json`, `phenotype_report.json` |
| `07_triage/` | `triage_decision.json`, `tier0_sourmash.json`, `local_triage.json` |
| `10_phylogeny/` | `tree.nwk`, `tree.png` |

---

## Test with public SRA data

```bash
# Install SRA tools: conda install -c bioconda sra-tools
fasterq-dump SRR32625477 --outdir data/raw_reads/ --threads 4
gzip data/raw_reads/SRR32625477_1.fastq data/raw_reads/SRR32625477_2.fastq

bash scripts/run_pipeline.sh --config workflow/test_config.yaml --cores 8
```

**Expected results for SRR32625477** (Haiti 2025, real WGS, ~7M reads):
- VRS: **33 / 🟢 LOW**
- Serotype: **Ogawa**, Toxin: **ctxB7**
- CTXφ: **NOT_DETECTED**, AMR genes: **0**
- Evo2 trajectory: **STABLE_ENDEMIC**
- Runtime: ~20 min on 8 cores / 16GB RAM

---

## Troubleshooting

| Issue | Fix |
|-------|-----|
| `WildcardError: NVIDIA_API_KEY` | Never put `"${VAR}"` in YAML — use `""` and pass via `--config` at runtime |
| Workspace locked (`.snakemake/locks/`) | `rm -f .snakemake/locks/*.lock` — stale lock from a killed run |
| OOM during Kraken2 | Use `kraken_db: "data/kraken2_serogroup"` or set `memory_management.force_tier: "BALANCED"` |
| Pipeline stops mid-run | Rerun — `run_pipeline.sh` always passes `--rerun-incomplete` |
| HyenaDNA not loading | See [docs/HYENADNA_SETUP.md](docs/HYENADNA_SETUP.md) |
| Low coverage warnings | Expected for <50k read samples — pipeline completes, just flags LOW_COVERAGE |

---

## Citation

If you use Vibrion Sentinel in published work:

> Vibrion Sentinel v2.0 — Genomic Surveillance for *Vibrio cholerae*.  
> https://github.com/intelogroup/vibrion-sentinel

---

## License

MIT License. See `LICENSE` for details.
