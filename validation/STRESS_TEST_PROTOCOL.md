# Vibrion Sentinel Field Stress Test Protocol
**Objective:** Validate computational predictions against biological reality before field deployment.

## 1. Hypervirulence & Toxin Profiler (rtxA Null Mutation Detection)

### Test Objective
Ensure the pipeline does not over- or under-estimate virulence by detecting the **G13602A stop codon in rtxA** (MARTX inactivation) specific to Haiti 2010 epidemic strains.

### Test Matrix
| Sample ID | Source | Expected rtxA Status | Expected Classification | Pass Criteria |
|-----------|--------|----------------------|------------------------|---|
| 2010EL-1786 | Haiti Reference (2010) | G13602A STOP (Inactive) | Hypervirulent (Hemolysin-dependent) | Correctly flags MARTX inactivation |
| HE-09 | Environmental variant | Functional (WT sequence) | Standard virulence | Identifies intact MARTX |
| EnvJ515 | Environmental (2018) | G13602A STOP (Inactive) | Haiti-like virulence | Matches ancestral profile |

### Implementation
- **Module:** `backend/core/logic/virulence_profiler.py`
- **Input:** `04_phylogeny/snp_calls.vcf` (or aligned consensus FASTA)
- **Output:** `02_serogroup/virulence_report.json` (new field: `rtxA_status`, `rtxA_mutation`, `rtxA_codon_position`)
- **Check Logic:**
  ```
  IF rtxA codon position 13602 is GAA → TAA (stop):
    CLASSIFY as "Haiti-like Hypervirulent (Hemolysin-dominant)"
  ELSE IF rtxA is intact:
    CLASSIFY as "Standard El Tor (MARTX-dependent)"
  ELSE:
    FLAG as "Unexpected rtxA deletion/rearrangement"
  ```

---

## 2. Time Capsule Clock Calibration (EnvJ515 Bridging)

### Test Objective
Confirm SNP distance metric correctly places environmental isolates as evolutionary bridges between 2010 ancestor and 2022 resurgence.

### Test Matrix
| Sample ID | Source | Collection Date | Expected SNP Distance to 2010 Ancestor | Expected Phylo Position | Pass Criteria |
|-----------|--------|-----------------|----------------------------------------|-------------------------|---|
| 2010EL-1786 | Haiti Reference | 2010 | 0 SNPs | Root/Ancestor | Returns 0 |
| EnvJ515 | Environmental | 2018 | 12-18 SNPs (estimated) | Basal node between 2010 & 2022 | Correctly places as bridge |
| Env5156 | Environmental | 2016 | 8-15 SNPs (estimated) | Intermediate | Temporal ordering respected |
| Env4303 | Environmental | 2015 | 10-16 SNPs (estimated) | Intermediate | Temporal ordering respected |
| 2022 Resurgent (example) | Clinical | 2022 | 27-29 SNPs | Derived clade | Matches literature distance |

### Implementation
- **Module:** `workflow/rules/phylogeny.smk` (IQ-TREE2 + TreeTime)
- **Input:** `04_phylogeny/aligned_msa.fasta`, `04_phylogeny/snp_calls.vcf`
- **Output:** `04_phylogeny/distance_metrics.json` (new fields: `snp_distance_to_reference`, `phylo_position`, `estimated_divergence_date`)
- **Check Logic:**
  ```
  SNP distance = count_variant_positions(sample_vs_ref)
  IF 2015-2018 samples cluster between 2010 & 2022:
    PASS (bridging confirmed)
  ELSE IF 2015-2018 samples appear as independent introductions:
    FAIL (clock uncalibrated)
  ELSE IF 2015-2018 samples group with 2022 (missing intermediate):
    PARTIAL (resolution insufficient)
  ```

---

## 3. Serology System (Inaba/Ogawa Switch via wbeT Frameshift)

### Test Objective
Detect loss-of-function mutations in *wbeT* gene (GAA→TAA stop codon) causing Ogawa→Inaba serotype transition.

### Test Matrix
| Sample ID | Source | Expected Serotype | wbeT Status | Expected Alert | Pass Criteria |
|-----------|--------|-------------------|-------------|---|---|
| 2010EL-1786 | Haiti Reference | Ogawa | WT (functional) | None | Correctly identifies WT |
| 2012EL-1410 | Clinical Inaba (Haiti) | Inaba | GAA→TAA STOP | "Vaccine Mismatch Risk" | Flags stop codon + serotype shift |
| HE-09 | Environmental Ogawa | Ogawa | WT | None | Identifies functional gene |

### Implementation
- **Module:** `backend/core/logic/serogroup_inference.py` (extend existing wbeT logic)
- **Input:** `04_phylogeny/snp_calls.vcf` (or k-mer evidence from `workflow/data/reference_kmers/wbeT_mutations.txt`)
- **Output:** `02_serogroup/serogroup_report.json` (extend field: `wbeT_mutation`, `wbeT_frameshift_detected`, `serotype_shift_alert`)
- **Check Logic:**
  ```
  IF wbeT codon region contains stop mutation (GAA→TAA at position X):
    SEROTYPE = Inaba
    ALERT = "Loss-of-function wbeT: Vaccine escape risk (serotype mismatch)"
  ELSE:
    SEROTYPE = Ogawa (or check rfb_c for Hikojima)
    ALERT = None
  ```

---

## 4. HGT & Plasmid Awareness (IncA/C vs SXT/ICE)

### Test Objective
Distinguish between chromosomal SXT/ICE resistance and plasmid-mediated IncA/C resistance (different transmission dynamics & regulators).

### Test Matrix
| Sample ID | Source | Primary AMR Element | Secondary Elements | Expected Discrimination | Pass Criteria |
|-----------|--------|---------------------|-------------------|---|---|
| 2010EL-1786 | Haiti Reference | SXT/R391 (chromosome) | None | "SXT-mediated (AcaCD regulation)" | Correctly identifies SXT |
| HC-36A1-like | Clinical (hypothetical) | IncA/C plasmid | Lacks SXT | "Plasmid-mediated (SetCD regulation)" | Distinguishes plasmid vs ICE |
| EnvJ515 | Environmental | SXT/R391 (chromosome) | None | "SXT-mediated" | Identifies chromosome-integrated |

### Implementation
- **Module:** `backend/core/logic/amr_profiler.py` (extend existing RGI + SXT logic)
- **Input:** `06_amr/rgi_report.json`, `05_structural/sxt_report.json`, assembly graph or plasmid predictions
- **Output:** `06_amr/amr_report.json` (new fields: `replicon_type`, `amr_element_location`, `regulator_set`)
- **Check Logic:**
  ```
  IF SXT element detected AND chromosome-integrated:
    ELEMENT_TYPE = "SXT/R391 (Integrative Conjugative Element)"
    REGULATORS = ["AcaCD", "SetCD" (if present)]
    RESISTANCE_GENES = [floR, tetR(D), trimethoprim, chloramphenicol] (SXT profile)
  ELSE IF Plasmid replicon detected (IncA/C signature):
    ELEMENT_TYPE = "IncA/C Conjugative Plasmid"
    REGULATORS = ["SetCD", "TraC" (if present)]
    RESISTANCE_GENES = [aac(3)-IIa, dfrA, tetA(D), floR] (plasmid profile)
    TRANSMISSION_DYNAMICS = "Higher horizontal gene transfer rate"
  ELSE:
    ELEMENT_TYPE = "Unclassified or mono-resistance"
  ```

---

## 5. Environmental Resilience (hapR & vpsA Integrity)

### Test Objective
Verify that *vpsA* and *hapR* are not only present but functionally intact (no frameshift/stop codons); biofilm phenotype correlates with regulatory integrity.

### Test Matrix
| Sample ID | Source | hapR Status | vpsA Status | Expected Biofilm Phenotype | Pass Criteria |
|-----------|--------|-------------|------------|---|---|
| 2010EL-1786 | Haiti Reference | WT (HapR competent) | WT (functional) | Rugose (37°C favored) | Flags both genes as WT |
| Hypothetical hapR mutant | Experimental | Stop codon | WT | Smooth (Loss of quorum sensing) | Detects hapR LoF |
| Hypothetical vpsA mutant | Experimental | WT | Frameshift | Smooth (Loss of capsule) | Detects vpsA LoF |

### Implementation
- **Module:** `backend/core/logic/environmental_profile.py` (extend existing vps cluster logic)
- **Input:** `04_phylogeny/snp_calls.vcf`, VCF-region typing for hapR/vpsA
- **Output:** `02_serogroup/serogroup_report.json` (new fields: `hapR_integrity`, `vpsA_integrity`, `biofilm_phenotype`)
- **Check Logic:**
  ```
  hapR_status = check_for_frameshift_or_stop(hapR_sequence, reference)
  vpsA_status = check_for_frameshift_or_stop(vpsA_sequence, reference)
  
  IF hapR_status == "WT" AND vpsA_status == "WT":
    BIOFILM_PHENOTYPE = "Rugose (HapR+, VpsA+)"
    RESILIENCE = "High (quorum sensing + matrix intact)"
  ELSE IF hapR_status == "LoF" OR vpsA_status == "LoF":
    BIOFILM_PHENOTYPE = "Smooth or reduced"
    RESILIENCE = "Reduced (key regulator/structural gene compromised)"
  ELSE:
    RESILIENCE = "Uncertain (partial mutations)"
  ```

---

## 6. Specificity / Imposter Detection (Bengal L1 vs L2 & Philippines GI-119)

### Test Objective
Reject foreign relatives (Bengal Lineage 1, Philippines "Hybrid El Tor") as non-endemic; distinguish via lineage-specific genomic islands.

### Test Matrix
| Sample ID | Source | Lineage | GI-119 Presence | Expected Rejection | Pass Criteria |
|-----------|--------|---------|---|---|---|
| 2010EL-1786 | Haiti Ancestral | Haiti-Yemen-Global | NO | Accept as endemic | Identifies as endemic |
| Bengal_L2_sample | Bangladesh (Lineage 2) | Bengal L2 (Global-linked) | NO | Accept (related to Haiti) | Allows L2 (Haiti-linked) |
| Bengal_L1_sample | Bangladesh (Lineage 1) | Bengal L1 (Dhaka-endemic) | NO (has unique R-M system) | REJECT as foreign lineage | Flags as non-endemic |
| Philippines_sample | Philippines outbreak | Hybrid El Tor | YES (GI-119) | REJECT as foreign lineage | Identifies GI-119 signature |

### Implementation
- **Module:** `backend/core/logic/regional_selection.py` & `global_matching.py`
- **Input:** `04_phylogeny/aligned_msa.fasta`, global k-mer reference database
- **Output:** `04_phylogeny/global_match.json` (new fields: `lineage_classification`, `foreign_lineage_flags`, `gi119_detected`, `restriction_modification_system`)
- **Check Logic:**
  ```
  IF sample clusters with Haiti-Yemen-Global in phylogeny:
    LINEAGE = "Haiti-Yemen-Global (L2)"
    STATUS = "Accept (endemic-linked)"
  ELSE IF sample clusters with Bengal L1 + lacks GI-119 + matches Dhaka-specific R-M:
    LINEAGE = "Bengal Lineage 1 (non-endemic)"
    STATUS = "REJECT (foreign lineage)"
    ALERT = "This is a local Dhaka strain, not the global pandemic lineage"
  ELSE IF sample contains GI-119 + IncA/C plasmid profile + Philippines k-mer signature:
    LINEAGE = "Philippines Hybrid El Tor (non-endemic)"
    STATUS = "REJECT (foreign lineage)"
    ALERT = "This is a distinct Philippines outbreak strain"
  ELSE:
    LINEAGE = "Unknown"
    STATUS = "Manual review required"
  ```

---

## 7. Sample Quality & Degradation (Freeze–Thaw Proxy)

### Test Objective
Define absolute limits of metagenomic triage under degraded DNA (cold chain failure, freeze–thaw cycles).

### Test Matrix
| Sample ID | Condition | Total Reads | Vibrio % | Mean k-mer CV | Expected QC Status | Pass Criteria |
|-----------|-----------|-------------|---------|---|---|---|
| Fresh extract | Immediate processing | 1M | 98.5% | 0.05 (low) | PASS | Purity >95%, CV <0.10 |
| Freeze-thaw 1x | 1 F-T cycle | 950k | 95.2% | 0.08 (moderate) | PASS | Purity >90%, CV <0.15 |
| Freeze-thaw 3x | 3 F-T cycles | 800k | 88.3% | 0.18 (high) | BORDERLINE | Purity 80-90%, CV >0.15 |
| Freeze-thaw 5x | 5 F-T cycles | 500k | 75.1% | 0.35 (very high) | FAIL | Purity <80%, CV >0.25 |

### Implementation
- **Module:** `backend/core/logic/qc_gating.py` (extend coverage integrity checks)
- **Input:** `00_qc/fastp.json`, `01_taxonomy/vibrio_stats.json`, k-mer coverage profile
- **Output:** `07_validation/checksum.json` (new fields: `degradation_proxy_cv`, `freeze_thaw_risk_score`, `dna_integrity_estimate`)
- **Check Logic:**
  ```
  cv_kmers = calculate_coefficient_of_variation(kmer_depths)
  degradation_score = (1 - vibrio_pct/100) + (cv_kmers / 0.5)  # normalized 0-2
  
  IF vibrio_pct > 95 AND cv_kmers < 0.10:
    QC_STATUS = "PASS (excellent quality, pristine sample)"
  ELSE IF vibrio_pct > 90 AND cv_kmers < 0.15:
    QC_STATUS = "PASS (acceptable quality, minor degradation)"
  ELSE IF vibrio_pct > 80 AND cv_kmers < 0.25:
    QC_STATUS = "BORDERLINE (low quality, significant degradation)"
    WARNING = "Freeze-thaw suspected; confidence in minor variants reduced"
  ELSE:
    QC_STATUS = "FAIL (excessive degradation, unreliable consensus)"
    ACTION = "Do not use for clinical decision-making"
  ```

---

## Stress Test Execution Plan

### Phase 1: Data Preparation
1. **Acquire or simulate test genomes:**
   - Haiti 2010EL-1786 (reference) — download from NCBI
   - HE-09 (environmental) — simulate from k-mer signature if unavailable
   - EnvJ515, Env5156, Env4303 (environmental bridges) — download or model SNP distances
   - 2012EL-1410 (Inaba) — download or construct synthetic variant
   - Bengal L1 & L2 samples — download from NCBI (BioSample metadata)
   - Philippines sample — download if available, else model GI-119 signature

2. **Organize in `data/validation_samples/`:**
   ```
   data/validation_samples/
   ├── haiti_2010el1786/          # Reference
   ├── env_he09/                  # Environmental variant
   ├── env_j515_2018/             # Bridge isolate
   ├── env_5156_2016/             # Bridge isolate
   ├── env_4303_2015/             # Bridge isolate
   ├── clinical_inaba_2012el1410/ # Serotype switch
   ├── bangladesh_l1/             # Foreign lineage
   ├── philippines_outbreak/      # Foreign lineage
   ├── fresh_extract/             # QC test (fresh)
   ├── freeze_thaw_1x/            # QC test (1 cycle)
   ├── freeze_thaw_3x/            # QC test (3 cycles)
   └── freeze_thaw_5x/            # QC test (5 cycles)
   ```

### Phase 2: Pipeline Execution
For each test sample:
```bash
snakemake -s workflow/Snakefile \
  --configfile workflow/validation_config.yaml \
  --config sample_id=<TEST_SAMPLE> \
  -j 4 \
  data/pipeline_output/<TEST_SAMPLE>/08_comprehensive_report/surveillance_report.md
```

### Phase 3: Output Parsing & Scoring
Use `validation/validation_harness.py` to:
- Extract JSON artifacts from each run
- Parse expected vs observed results for all 7 objectives
- Assign PASS/FAIL/PARTIAL per test
- Generate `STRESS_TEST_RESULTS.json` with detailed metrics

### Phase 4: Gap Documentation
List which objectives require additional implementation or have edge-case failures.

---

## Success Criteria

| Objective | PASS Threshold | Current Status | Gap Notes |
|-----------|---|---|---|
| 1. rtxA Stop-Codon Detection | Correctly identifies G13602A in Haiti ref; identifies WT in HE-09 | Unknown (not implemented) | Requires addition to virulence profiler |
| 2. Time Capsule Clock | Correct SNP distances for 2015–2018 bridges; proper phylo ordering | Partial (TreeTime runs but EnvJ515 not tested) | Needs explicit temporal validation dataset |
| 3. wbeT Frameshift Detection | Flags Inaba stop codon; alerts on vaccine mismatch | Partial (wbeT presence detected, but not frameshift) | Requires frameshift/stop-codon parser |
| 4. IncA/C vs SXT Discrimination | Distinguishes plasmid from SXT; identifies regulators | Minimal (SXT only; no plasmid logic) | Requires RGI + replicon integration |
| 5. hapR/vpsA Integrity | Flags LoF mutations; correlates with biofilm phenotype | Minimal (genes present in signatures, not integrity checked) | Requires LoF detection logic |
| 6. Lineage Specificity | Rejects Bengal L1 & Philippines; accepts L2 | Minimal (global matching exists; discriminators not specific) | Requires GI-119 & lineage-specific classifiers |
| 7. Degradation Proxy | QC gates fail at CV >0.25 & vibrio <80%; passes fresh samples | Partial (purity checked, CV proxy not implemented) | Requires k-mer CV metric + freeze-thaw risk scoring |

---

## Outputs

On completion, generate:
1. **`STRESS_TEST_RESULTS.json`** — Detailed per-sample metrics for all 7 objectives
2. **`STRESS_TEST_RESULTS.md`** — Human-readable summary with recommendations
3. **Field Deployment Confidence Statement** — Can field teams use this in production? What caveats apply?

