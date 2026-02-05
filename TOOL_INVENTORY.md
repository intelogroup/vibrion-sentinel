# Vibrion Sentinel Pipeline - Complete Tool Inventory

## Summary
After analyzing the vibrion-public pipeline, we have identified **30+ specialized bioinformatics tools** organized into 4 functional categories. All tools are now documented in the README.md architecture section.

---

## Pipeline Stages & Tools

### Stage 0: Quality Control
- **Fastp** (v0.23.0+): Quality filtering with sliding window (Q20, 4bp window)

### Stage 1: Decontamination  
- **Hostile** (v2.0.2+): Aggressive host DNA removal
- **Minimap2** (v2.26+): Short-read mapping for filtering

### Stage 2: Species Classification
- **Kraken2** (v2.1.3+): k-mer based classification
  - Custom Haiti database
  - Standard 8GB database
  - Viral database (backup)
  - Phage sentinel module

### Stage 3: Read Rescue
- **MMseqs2**: Protein-based alignment for mutated reads

### Stage 4: Alignment
- **BWA** (v0.7.17+): Burrows-Wheeler aligner for reference mapping
- **Minimap2** (v2.26+): Long/short-read aligner

### Stage 5: De Novo Assembly (Fallback)
- **SPAdes** (v3.15+): Assembly when coverage insufficient

### Stage 6: Consensus Generation
- **Samtools** (v1.18+): SAM/BAM processing and pileup
- **Pilon** (v1.24+): Illumina-based consensus refinement

### Stage 7: Genome Polishing
**Nanopore path:**
- **Medaka** (v1.11.3): Nanopore base calling and error correction
- **Minimap2**: Alignment for nanopore reads
  
**Illumina path:**
- **Pilon** (v1.24): Illumina-based polishing
- **BWA** (v0.7.17): Short-read alignment

**Cross-platform correction:**
- **Polypolish** (v0.6.0): Multiple alignment polishing
- **FMLRC2** (v0.1.7): Long-read error correction

### Stage 8: Multiple Sequence Alignment
- **MAFFT** (v7.50+): Fast phylogenetic alignment

### Stage 9: Variant Calling
- **BCFtools** (v1.18+): VCF calling from BAM files
- **Freebayes** (v1.3.7+): Bayesian SNP/indel detection

### Stage 10: Annotation
- **SnpEff** (v5.1+): Variant functional annotation
- **BLAST** (v2.14.1+): Sequence similarity and AMR detection

---

## Intelligence Tier Tools

### Tier 0: Flash Triage (k-mer based)
- **Sourmash** (v4.8.0+): MinHash sketching and sequence comparison
- **Screed** (v1.0.0+): Fast sequence metadata storage

### Tier 1: Local AI Analysis
- **HyenaDNA**: Structural anomaly detection (local CPU model)
- **PyTorch** (v2.0.0+): Deep learning runtime
- **Transformers** (v4.35.0+): Pre-trained models

### Tier 2: Cloud Deep Forensics
- **Evo2** (cloud-based): Escalation analysis

---

## Specialized Detection Tools

### AMR Detection
- **RGI** (Resistance Gene Identifier): AMR phenotype prediction
- **DIAMOND** (v0.9.36+): Fast protein alignment
- **Prodigal**: Gene prediction for ORF detection
- **BLAST** (v2.14.1+): Sequence comparison

### Phylogenetics & Evolution
- **IQ-TREE** (v2.2.0+): Fast phylogenetic inference
- **TreeTime** (v0.11.1+): Molecular clock analysis
- **Augur**: Bioconductor phylogenetic toolkit
- **Bioconductor** (R packages):
  - ggtree: Tree visualization
  - treeio: Tree I/O
  - ape: Analysis of phylogenetics
  - ggplot2: Publication-quality plots

### Data Transport & Async
- **Aiohttp** (v3.9+): Async HTTP for API calls
- **Httpx** (v0.27+): Modern HTTP client
- **Requests** (v2.31+): HTTP library

### Utilities
- **Pigz** (v2.8+): Parallel gzip compression
- **Dnaio**: DNA I/O library (hostile dependency)
- **Biopython** (v1.80+): Sequence parsing and manipulation
- **Pysam** (v0.20+): SAM/BAM interface
- **NumPy** (v1.24+): Numerical computing
- **Matplotlib** (v3.8.0+): Plotting and visualization

---

## Environment Organization

The pipeline uses **8 separate conda environments** for modularity:

1. **analysis.yaml** - Main analysis tools (fastp, mafft, bwa, samtools, bcftools, PyTorch, transformers)
2. **assembly.yaml** - SPAdes, samtools, biopython
3. **polishing_illumina.yaml** - Pilon, BWA, samtools, bcftools (Illumina-optimized)
4. **polishing_nanopore.yaml** - Medaka, minimap2, samtools, bcftools (Nanopore-optimized)
5. **phylogeny.yaml** - IQ-TREE, MAFFT, TreeTime, Augur, R packages
6. **kraken2.yaml** - Kraken2, pigz
7. **hostile.yaml** - Hostile, minimap2
8. **amr_rgi.yaml** - RGI, DIAMOND, Prodigal, BLAST

---

## Updates Made to README.md

✅ Added comprehensive 3-part architecture section:
- **Part 1**: 11-stage data pipeline with all tools
- **Part 2**: 3-tier intelligence system
- **Part 3**: Supporting tools table

✅ Included version requirements where available
✅ Clarified platform-specific tooling (Nanopore vs Illumina)
✅ Added acknowledgment to CholeraSeq pipeline with official documentation link

---

## Recommended Next Steps

1. **Version pinning**: Review minimum version requirements for each tool
2. **Docker optimization**: Ensure all tools properly containerized in Dockerfile
3. **Performance profiling**: Document tool runtime and memory requirements per stage
4. **Alternative tools**: Document fallback tools for each critical stage
5. **Validation tests**: Add sample runs with expected tool outputs

