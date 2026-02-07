# Snakemake Pipeline Integration - Complete

**Date:** February 7, 2026  
**Status:** ✅ STRESS TEST MODULES SUCCESSFULLY INTEGRATED INTO PIPELINE

---

## Summary

All 7 biological validation modules have been successfully wired into the Snakemake pipeline. The stress test validation framework is now operational and will execute automatically during pipeline runs.

---

## Integration Details

### Files Modified

1. **`workflow/rules/stress_test.smk`** (NEW)
   - Comprehensive rule: `stress_test_all`
   - Orchestrates all 7 biological validation modules
   - Generates 8 output files with validation results

2. **`workflow/Snakefile`**
   - Added include statement: `include: "rules/stress_test.smk"`
   - Added stress test outputs to `rule all` targets

### Rule: `stress_test_all`

**Input Requirements:**
- Annotated VCF (from SNP calling + SnpEff annotation)
- BAM file (from alignment to reference)
- Reference FASTA sequence

**Outputs Generated:**

1. `02_serogroup/virulence_report.json` — rtxA G13602A detection
2. `02_serogroup/serotype_mutations.json` — wbeT frameshift + vaccine mismatch
3. `06_amr/amr_element_classification.json` — IncA/C vs SXT/ICE discrimination
4. `05_resilience/biofilm_phenotype.json` — hapR/vpsA integrity + phenotype
5. `04_phylogeny/lineage_classification.json` — Lineage classification (Haiti-L2, etc.)
6. `07_validation/degradation_metrics.json` — k-mer CV, SNP distance, freeze-thaw estimate
7. `07_validation/stress_test_comprehensive.json` — Combined results from all 6 modules
8. `07_validation/stress_test_summary.md` — Human-readable Markdown report

---

## Pipeline Execution

### Run Full Pipeline with Stress Tests

```bash
# All samples configured in config.yaml
snakemake --configfile workflow/config/config.yaml --cores 4

# Single sample only
snakemake --configfile workflow/config/config.yaml \
  --cores 4 \
  data/pipeline_output/sample_name/07_validation/stress_test_summary.md
```

### Expected Outputs

Each sample will have:
```
data/pipeline_output/{sample}/
  ├── 02_serogroup/
  │   ├── virulence_report.json
  │   └── serotype_mutations.json
  ├── 04_phylogeny/
  │   └── lineage_classification.json
  ├── 05_resilience/
  │   └── biofilm_phenotype.json
  ├── 06_amr/
  │   └── amr_element_classification.json
  └── 07_validation/
      ├── degradation_metrics.json
      ├── stress_test_comprehensive.json
      └── stress_test_summary.md
```

---

## Field Deployment Validation

### Deployment Decision Logic

Field deployment ready when ALL conditions met:

```python
field_deployment_ready = (
    lineage_classification in ["Haiti-L2", "Yemen-L2", "Global-L2"] AND
    sample_quality_pass (kmer_cv < 0.25) AND
    NOT vaccine_mismatch_alert
)
```

### Review Checklist

Before deploying to field, verify:

- ✅ All 7 objective outputs present
- ✅ No errors in log files (`logs/stress_test_all.log`)
- ✅ `field_deployment_ready` = true in `stress_test_comprehensive.json`
- ✅ `stress_test_summary.md` shows "✅ PASS" recommendation
- ✅ Sample quality metrics pass thresholds

---

## Test Execution

To verify the integration works:

```bash
# Dry-run to check rule dependencies
snakemake --configfile workflow/config/config.yaml --dry-run | grep stress_test

# Actual run (requires full pipeline completion first)
snakemake --configfile workflow/config/config.yaml --cores 4
```

---

## Module Integration Points

### How Modules Are Called

Each module is invoked from the consolidated `stress_test_all` rule:

1. **VirulenceProfiler** — Scans reference and VCF for rtxA mutations
2. **SerotypeMutationDetector** — Parses VCF for wbeT frameshifts
3. **AMRElementDiscriminator** — Identifies plasmid vs ICE signatures
4. **EnvironmentalResilienceProfiler** — Checks hapR/vpsA LoF mutations
5. **LineageSpecificityClassifier** — Compares against lineage markers
6. **DegradationProxyCalculator** — Calculates k-mer CV from BAM
7. **SNPDistanceCalculator** — Counts SNP distance to reference

### Data Flow

```
VCF + BAM + Reference
    ↓
stress_test_all rule (orchestrator)
    ├→ Module 1 → virulence_report.json
    ├→ Module 2 → serotype_mutations.json
    ├→ Module 3 → amr_element_classification.json
    ├→ Module 4 → biofilm_phenotype.json
    ├→ Module 5 → lineage_classification.json
    ├→ Module 6 → degradation_metrics.json
    └→ Aggregator
        ├→ stress_test_comprehensive.json
        └→ stress_test_summary.md
```

---

## GitHub Status

**Repository:** https://github.com/intelogroup/vibrion-sentinel  
**Latest Commit:** `4ba2661` - "feat: Integrate stress test validation rules into Snakemake pipeline"  
**Branch:** main  
**Status:** ✅ Pushed and synced

---

## Next Steps: Field Test Execution

### Phase 1: Validation Sample Set (1-2 weeks)
- Acquire or generate 8-12 test samples representing:
  - Haiti 2010EL-1786 (reference, G13602A)
  - HE-09 (environmental, WT virulence)
  - EnvJ515 (2018 environmental bridge)
  - 2012EL-1410 (clinical Inaba)
  - Bengali L1 strain (should be REJECTED)
  - Philippines GI-119 strain (should be REJECTED)
  - Freeze-thaw samples (0x, 1x, 3x, 5x cycles)

### Phase 2: Pipeline Execution
- Run full pipeline: `snakemake --cores 4 --configfile config.yaml`
- Collect stress test outputs for all samples
- Run validation harness: `python3 validation/validation_harness.py`

### Phase 3: Results Analysis
- Review `STRESS_TEST_RESULTS.md` + `.json`
- Verify all objectives PASS for appropriate samples
- Document any issues for field team briefing

### Phase 4: Field Deployment
- Brief field teams on stress test results
- Provide deployment checklist
- Deploy to Haiti public health facilities

---

## Documentation References

- **STRESS_TEST_PROTOCOL.md** — Detailed test matrices and expectations
- **IMPLEMENTATION_GUIDE.md** — Module implementation details
- **README.md** — Quick start guide for field teams
- **VIBRION_PUBLIC_TEST_REPORT.md** — Current test status

---

## Support

For issues during field execution:

1. Check `logs/stress_test_all.log` for error messages
2. Verify module imports: `python3 backend/core/logic/{module}.py`
3. Ensure input files (VCF, BAM) exist and are valid
4. Consult `IMPLEMENTATION_GUIDE.md` for troubleshooting

---

✅ **Snakemake integration complete. Pipeline ready for field validation testing.**
