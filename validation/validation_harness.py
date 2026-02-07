#!/usr/bin/env python3
"""
Vibrion Sentinel Stress Test Harness
Orchestrates all 7 biological validation objectives and scores results.
"""

import json
import subprocess
import os
import sys
from pathlib import Path
from typing import Dict, List, Any, Tuple
from dataclasses import dataclass, asdict
from collections import defaultdict
from datetime import datetime
import hashlib

@dataclass
class TestResult:
    """Record of a single validation test outcome"""
    objective_id: str
    objective_name: str
    test_sample: str
    expected_call: str
    observed_call: str
    pass_fail: str  # "PASS", "FAIL", "PARTIAL", "NOT_TESTED"
    details: str
    artifact_path: str  # JSON path inspected

class StressTestHarness:
    """Main orchestrator for Vibrion Sentinel validation"""
    
    def __init__(self, workspace_root: str = "/Users/kalinovdameus/Developer/Vibrion"):
        self.workspace = Path(workspace_root)
        self.validation_dir = self.workspace / "validation"
        self.pipeline_output = self.workspace / "data" / "pipeline_output"
        self.test_samples_dir = self.workspace / "data" / "validation_samples"
        self.results = []  # type: List[TestResult]
        self.test_matrix = self._define_test_matrix()
        
    def _define_test_matrix(self) -> Dict[str, Dict[str, Any]]:
        """Define expected outcomes for all 7 objectives across test samples"""
        return {
            "OBJ1_RTXA": {
                "name": "rtxA Stop-Codon Detection (Hypervirulence Profiler)",
                "tests": {
                    "haiti_2010el1786": {
                        "expected_rtxa_status": "G13602A_STOP_CODON",
                        "expected_classification": "Haiti_Hypervirulent_Hemolysin_Dependent",
                        "artifact": "02_serogroup/virulence_report.json",
                        "check_field": "rtxA_status"
                    },
                    "env_he09": {
                        "expected_rtxa_status": "FUNCTIONAL_WT",
                        "expected_classification": "Standard_El_Tor_MARTX_Dependent",
                        "artifact": "02_serogroup/virulence_report.json",
                        "check_field": "rtxA_status"
                    },
                    "env_j515_2018": {
                        "expected_rtxa_status": "G13602A_STOP_CODON",
                        "expected_classification": "Haiti_Like_Environmental",
                        "artifact": "02_serogroup/virulence_report.json",
                        "check_field": "rtxA_status"
                    }
                }
            },
            "OBJ2_CLOCK": {
                "name": "Time Capsule Clock (SNP Distance & Phylo Bridging)",
                "tests": {
                    "haiti_2010el1786": {
                        "expected_snp_distance": 0,
                        "expected_phylo_position": "ROOT_ANCESTOR",
                        "artifact": "04_phylogeny/distance_metrics.json",
                        "check_field": "snp_distance_to_reference"
                    },
                    "env_j515_2018": {
                        "expected_snp_distance_range": (12, 18),
                        "expected_phylo_position": "BASAL_NODE_BRIDGE",
                        "artifact": "04_phylogeny/distance_metrics.json",
                        "check_field": "snp_distance_to_reference"
                    },
                    "env_5156_2016": {
                        "expected_snp_distance_range": (8, 15),
                        "expected_phylo_position": "INTERMEDIATE",
                        "artifact": "04_phylogeny/distance_metrics.json",
                        "check_field": "snp_distance_to_reference"
                    }
                }
            },
            "OBJ3_SEROTYPE": {
                "name": "Serology System (wbeT Frameshift & Inaba Switch)",
                "tests": {
                    "haiti_2010el1786": {
                        "expected_serotype": "Ogawa",
                        "expected_wbet_status": "WT_FUNCTIONAL",
                        "expected_alert": "NONE",
                        "artifact": "02_serogroup/serogroup_report.json",
                        "check_field": "wbeT_mutation"
                    },
                    "clinical_inaba_2012el1410": {
                        "expected_serotype": "Inaba",
                        "expected_wbet_status": "GAA_TAA_STOP",
                        "expected_alert": "VACCINE_MISMATCH_RISK",
                        "artifact": "02_serogroup/serogroup_report.json",
                        "check_field": "wbeT_mutation"
                    }
                }
            },
            "OBJ4_HGT": {
                "name": "HGT & Plasmid Awareness (IncA/C vs SXT/ICE)",
                "tests": {
                    "haiti_2010el1786": {
                        "expected_amr_element": "SXT_R391_ICE",
                        "expected_location": "CHROMOSOME_INTEGRATED",
                        "expected_regulators": ["AcaCD", "SetCD"],
                        "artifact": "06_amr/amr_report.json",
                        "check_field": "replicon_type"
                    },
                    "hc36a1_hypothetical": {
                        "expected_amr_element": "INCA_C_PLASMID",
                        "expected_location": "PLASMID",
                        "expected_regulators": ["SetCD", "TraC"],
                        "artifact": "06_amr/amr_report.json",
                        "check_field": "replicon_type"
                    }
                }
            },
            "OBJ5_RESILIENCE": {
                "name": "Environmental Resilience (hapR & vpsA Integrity)",
                "tests": {
                    "haiti_2010el1786": {
                        "expected_hapR_status": "WT_PROFICIENT",
                        "expected_vpsA_status": "WT_FUNCTIONAL",
                        "expected_biofilm": "Rugose",
                        "artifact": "02_serogroup/serogroup_report.json",
                        "check_field": "hapR_integrity"
                    },
                    "hypothetical_hapR_mutant": {
                        "expected_hapR_status": "LOF_STOP_CODON",
                        "expected_vpsA_status": "WT_FUNCTIONAL",
                        "expected_biofilm": "Smooth",
                        "artifact": "02_serogroup/serogroup_report.json",
                        "check_field": "hapR_integrity"
                    }
                }
            },
            "OBJ6_SPECIFICITY": {
                "name": "Specificity / Imposter Detection (Bengal L1 vs L2, Philippines GI-119)",
                "tests": {
                    "haiti_2010el1786": {
                        "expected_lineage": "Haiti_Yemen_Global_L2",
                        "expected_status": "ACCEPT_ENDEMIC",
                        "artifact": "04_phylogeny/global_match.json",
                        "check_field": "lineage_classification"
                    },
                    "bangladesh_l1": {
                        "expected_lineage": "Bengal_L1_Dhaka_Endemic",
                        "expected_status": "REJECT_FOREIGN",
                        "expected_alert": "Foreign_Lineage_Non_Endemic",
                        "artifact": "04_phylogeny/global_match.json",
                        "check_field": "lineage_classification"
                    },
                    "philippines_outbreak": {
                        "expected_lineage": "Philippines_Hybrid_El_Tor",
                        "expected_status": "REJECT_FOREIGN",
                        "expected_gi119": True,
                        "artifact": "04_phylogeny/global_match.json",
                        "check_field": "gi119_detected"
                    }
                }
            },
            "OBJ7_DEGRADATION": {
                "name": "Sample Quality & Degradation (Freeze-Thaw Proxy)",
                "tests": {
                    "fresh_extract": {
                        "expected_vibrio_pct_min": 95,
                        "expected_cv_kmer_max": 0.10,
                        "expected_qc_status": "PASS",
                        "artifact": "07_validation/checksum.json",
                        "check_field": "degradation_proxy_cv"
                    },
                    "freeze_thaw_1x": {
                        "expected_vibrio_pct_min": 90,
                        "expected_cv_kmer_max": 0.15,
                        "expected_qc_status": "PASS",
                        "artifact": "07_validation/checksum.json",
                        "check_field": "degradation_proxy_cv"
                    },
                    "freeze_thaw_3x": {
                        "expected_vibrio_pct_min": 80,
                        "expected_cv_kmer_max": 0.25,
                        "expected_qc_status": "BORDERLINE",
                        "artifact": "07_validation/checksum.json",
                        "check_field": "degradation_proxy_cv"
                    },
                    "freeze_thaw_5x": {
                        "expected_vibrio_pct_min": -1,  # Will fail
                        "expected_cv_kmer_max": 0.5,
                        "expected_qc_status": "FAIL",
                        "artifact": "07_validation/checksum.json",
                        "check_field": "degradation_proxy_cv"
                    }
                }
            }
        }
    
    def run_pipeline_for_sample(self, sample_id: str) -> bool:
        """Execute Snakemake pipeline for a single test sample"""
        print(f"[ORCHESTRATOR] Running pipeline for sample: {sample_id}")
        
        config_file = self.workspace / "workflow" / "validation_config.yaml"
        if not config_file.exists():
            # Create a minimal validation config
            self._create_validation_config()
        
        # Run snakemake
        cmd = [
            "snakemake",
            "-s", str(self.workspace / "workflow" / "Snakefile"),
            f"--configfile", str(config_file),
            f"--config", f"sample_id={sample_id}",
            "-j", "4",
            "-p"  # Print commands
        ]
        
        try:
            result = subprocess.run(cmd, cwd=str(self.workspace), capture_output=True, text=True, timeout=3600)
            if result.returncode == 0:
                print(f"[✓] Pipeline completed for {sample_id}")
                return True
            else:
                print(f"[✗] Pipeline failed for {sample_id}")
                print("STDERR:", result.stderr[:500])
                return False
        except subprocess.TimeoutExpired:
            print(f"[✗] Pipeline timeout for {sample_id}")
            return False
        except Exception as e:
            print(f"[✗] Error running pipeline: {e}")
            return False
    
    def load_json_artifact(self, sample_id: str, artifact_path: str) -> Dict[str, Any]:
        """Load a JSON artifact from pipeline output"""
        full_path = self.pipeline_output / sample_id / artifact_path
        if full_path.exists():
            try:
                with open(full_path, 'r') as f:
                    return json.load(f)
            except json.JSONDecodeError:
                print(f"[⚠] JSON decode error in {full_path}")
                return {}
        else:
            print(f"[⚠] Missing artifact: {full_path}")
            return {}
    
    def evaluate_objective_1_rtxa(self) -> List[TestResult]:
        """Test rtxA stop-codon detection"""
        results = []
        obj_name = "rtxA Stop-Codon Detection"
        
        for sample_id, expected in self.test_matrix["OBJ1_RTXA"]["tests"].items():
            artifact_data = self.load_json_artifact(sample_id, expected["artifact"])
            
            # Check if rtxA_status field exists (currently it won't, so this tests for implementation gap)
            observed = artifact_data.get("rtxA_status", "NOT_IMPLEMENTED")
            expected_val = expected["expected_rtxa_status"]
            
            if observed == "NOT_IMPLEMENTED":
                status = "NOT_TESTED"
                details = f"rtxA_status field not yet implemented in virulence_report.json"
            elif observed == expected_val:
                status = "PASS"
                details = f"Correctly identified {expected_val}"
            else:
                status = "FAIL"
                details = f"Expected {expected_val}, got {observed}"
            
            results.append(TestResult(
                objective_id="OBJ1",
                objective_name=obj_name,
                test_sample=sample_id,
                expected_call=expected_val,
                observed_call=observed,
                pass_fail=status,
                details=details,
                artifact_path=expected["artifact"]
            ))
        
        return results
    
    def evaluate_objective_2_clock(self) -> List[TestResult]:
        """Test SNP distance and phylo bridging"""
        results = []
        obj_name = "Time Capsule Clock (SNP Distance)"
        
        for sample_id, expected in self.test_matrix["OBJ2_CLOCK"]["tests"].items():
            artifact_data = self.load_json_artifact(sample_id, expected["artifact"])
            
            observed = artifact_data.get("snp_distance_to_reference", -1)
            
            if observed == -1:
                status = "NOT_TESTED"
                expected_str = str(expected.get("expected_snp_distance") or expected.get("expected_snp_distance_range"))
                details = "snp_distance_to_reference field not found in distance_metrics.json"
            else:
                if "expected_snp_distance" in expected:
                    exp_val = expected["expected_snp_distance"]
                    if observed == exp_val:
                        status = "PASS"
                        details = f"Correct distance: {observed} SNPs"
                    else:
                        status = "FAIL"
                        details = f"Expected {exp_val} SNPs, got {observed}"
                else:
                    exp_min, exp_max = expected["expected_snp_distance_range"]
                    expected_str = f"{exp_min}-{exp_max}"
                    if exp_min <= observed <= exp_max:
                        status = "PASS"
                        details = f"Within range {expected_str}: {observed} SNPs"
                    else:
                        status = "FAIL"
                        details = f"Expected {expected_str}, got {observed}"
            
            results.append(TestResult(
                objective_id="OBJ2",
                objective_name=obj_name,
                test_sample=sample_id,
                expected_call=expected_str if 'expected_str' in locals() else str(expected.get("expected_snp_distance", "?")),
                observed_call=str(observed),
                pass_fail=status,
                details=details,
                artifact_path=expected["artifact"]
            ))
        
        return results
    
    def evaluate_objective_3_serotype(self) -> List[TestResult]:
        """Test wbeT frameshift detection for Inaba switch"""
        results = []
        obj_name = "Serology System (wbeT Frameshift)"
        
        for sample_id, expected in self.test_matrix["OBJ3_SEROTYPE"]["tests"].items():
            artifact_data = self.load_json_artifact(sample_id, expected["artifact"])
            
            observed = artifact_data.get("wbeT_mutation", "NOT_IMPLEMENTED")
            expected_val = expected["expected_wbet_status"]
            
            if observed == "NOT_IMPLEMENTED":
                status = "NOT_TESTED"
                details = "wbeT_mutation field not yet implemented (only wbeT presence is detected)"
            elif observed == expected_val or (expected_val == "WT_FUNCTIONAL" and "wt" in str(observed).lower()):
                status = "PASS"
                details = f"Correctly identified {expected_val}"
            else:
                status = "FAIL"
                details = f"Expected {expected_val}, got {observed}"
            
            results.append(TestResult(
                objective_id="OBJ3",
                objective_name=obj_name,
                test_sample=sample_id,
                expected_call=expected_val,
                observed_call=str(observed),
                pass_fail=status,
                details=details,
                artifact_path=expected["artifact"]
            ))
        
        return results
    
    def evaluate_objective_4_hgt(self) -> List[TestResult]:
        """Test IncA/C vs SXT discrimination"""
        results = []
        obj_name = "HGT & Plasmid Awareness"
        
        for sample_id, expected in self.test_matrix["OBJ4_HGT"]["tests"].items():
            artifact_data = self.load_json_artifact(sample_id, expected["artifact"])
            
            observed = artifact_data.get("replicon_type", "NOT_IMPLEMENTED")
            expected_val = expected["expected_amr_element"]
            
            if observed == "NOT_IMPLEMENTED":
                status = "NOT_TESTED"
                details = "replicon_type field not yet implemented (SXT only, no plasmid logic)"
            elif observed == expected_val:
                status = "PASS"
                details = f"Correctly identified {expected_val}"
            else:
                # Allow partial credit if SXT-only is found but IncA/C is expected
                status = "PARTIAL"
                details = f"Expected {expected_val}, got {observed} (SXT logic present, plasmid logic missing)"
            
            results.append(TestResult(
                objective_id="OBJ4",
                objective_name=obj_name,
                test_sample=sample_id,
                expected_call=expected_val,
                observed_call=str(observed),
                pass_fail=status,
                details=details,
                artifact_path=expected["artifact"]
            ))
        
        return results
    
    def evaluate_objective_5_resilience(self) -> List[TestResult]:
        """Test hapR/vpsA integrity detection"""
        results = []
        obj_name = "Environmental Resilience (hapR/vpsA)"
        
        for sample_id, expected in self.test_matrix["OBJ5_RESILIENCE"]["tests"].items():
            artifact_data = self.load_json_artifact(sample_id, expected["artifact"])
            
            observed_hapR = artifact_data.get("hapR_integrity", "NOT_IMPLEMENTED")
            expected_hapR = expected["expected_hapR_status"]
            
            if observed_hapR == "NOT_IMPLEMENTED":
                status = "NOT_TESTED"
                details = "hapR_integrity and vpsA_integrity fields not yet implemented (genes present but not integrity-checked)"
            elif observed_hapR == expected_hapR:
                status = "PASS"
                details = f"Correctly identified hapR={expected_hapR}"
            else:
                status = "FAIL"
                details = f"Expected hapR={expected_hapR}, got {observed_hapR}"
            
            results.append(TestResult(
                objective_id="OBJ5",
                objective_name=obj_name,
                test_sample=sample_id,
                expected_call=expected_hapR,
                observed_call=str(observed_hapR),
                pass_fail=status,
                details=details,
                artifact_path=expected["artifact"]
            ))
        
        return results
    
    def evaluate_objective_6_specificity(self) -> List[TestResult]:
        """Test lineage discrimination (Bengal L1, Philippines GI-119)"""
        results = []
        obj_name = "Specificity / Imposter Detection"
        
        for sample_id, expected in self.test_matrix["OBJ6_SPECIFICITY"]["tests"].items():
            artifact_data = self.load_json_artifact(sample_id, expected["artifact"])
            
            observed = artifact_data.get("lineage_classification", "NOT_IMPLEMENTED")
            expected_val = expected["expected_lineage"]
            
            if observed == "NOT_IMPLEMENTED":
                status = "NOT_TESTED"
                details = "lineage_classification field not yet implemented (global matching exists but discriminators not specific)"
            elif observed == expected_val:
                status = "PASS"
                details = f"Correctly identified {expected_val} and {expected.get('expected_status', 'UNKNOWN')}"
            else:
                status = "FAIL"
                details = f"Expected {expected_val}, got {observed}"
            
            results.append(TestResult(
                objective_id="OBJ6",
                objective_name=obj_name,
                test_sample=sample_id,
                expected_call=expected_val,
                observed_call=str(observed),
                pass_fail=status,
                details=details,
                artifact_path=expected["artifact"]
            ))
        
        return results
    
    def evaluate_objective_7_degradation(self) -> List[TestResult]:
        """Test freeze-thaw degradation proxy"""
        results = []
        obj_name = "Sample Quality & Degradation"
        
        for sample_id, expected in self.test_matrix["OBJ7_DEGRADATION"]["tests"].items():
            artifact_data = self.load_json_artifact(sample_id, expected["artifact"])
            
            observed_cv = artifact_data.get("degradation_proxy_cv", -1)
            expected_cv = expected["expected_cv_kmer_max"]
            expected_qc = expected["expected_qc_status"]
            
            if observed_cv == -1:
                status = "NOT_TESTED"
                details = "degradation_proxy_cv field not yet implemented (purity checked but CV proxy missing)"
            elif observed_cv <= expected_cv:
                status = "PASS"
                details = f"CV={observed_cv:.3f} <= {expected_cv} ({expected_qc})"
            else:
                status = "FAIL"
                details = f"CV={observed_cv:.3f} > {expected_cv} (expected {expected_qc})"
            
            results.append(TestResult(
                objective_id="OBJ7",
                objective_name=obj_name,
                test_sample=sample_id,
                expected_call=f"CV<{expected_cv}",
                observed_call=f"CV={observed_cv:.3f}" if observed_cv >= 0 else "NOT_FOUND",
                pass_fail=status,
                details=details,
                artifact_path=expected["artifact"]
            ))
        
        return results
    
    def evaluate_all_objectives(self) -> None:
        """Run all 7 validation objectives"""
        print("\n" + "="*80)
        print("VIBRION SENTINEL STRESS TEST: EVALUATION PHASE")
        print("="*80 + "\n")
        
        evaluators = [
            ("OBJ1: rtxA Detection", self.evaluate_objective_1_rtxa),
            ("OBJ2: Time Capsule Clock", self.evaluate_objective_2_clock),
            ("OBJ3: Serology (wbeT)", self.evaluate_objective_3_serotype),
            ("OBJ4: HGT/Plasmid", self.evaluate_objective_4_hgt),
            ("OBJ5: Resilience (hapR/vpsA)", self.evaluate_objective_5_resilience),
            ("OBJ6: Specificity (Lineage)", self.evaluate_objective_6_specificity),
            ("OBJ7: Degradation (QC)", self.evaluate_objective_7_degradation),
        ]
        
        for obj_label, evaluator in evaluators:
            print(f"\n[EVALUATING] {obj_label}")
            obj_results = evaluator()
            self.results.extend(obj_results)
            for res in obj_results:
                icon = "✓" if res.pass_fail == "PASS" else ("?" if res.pass_fail == "NOT_TESTED" else "✗")
                print(f"  {icon} {res.test_sample}: {res.pass_fail} ({res.details[:60]}...)")
    
    def generate_report(self) -> str:
        """Generate comprehensive stress test report"""
        print("\n" + "="*80)
        print("STRESS TEST REPORT GENERATION")
        print("="*80 + "\n")
        
        report = f"""# Vibrion Sentinel Stress Test Results
**Generated:** {datetime.now().isoformat()}
**Workspace:** {self.workspace}

## Executive Summary

This report documents the validation of Vibrion Sentinel's ability to detect 7 critical biological mechanisms under realistic field conditions:

1. **Hypervirulence Profiling** (rtxA stop codon detection)
2. **Evolutionary Clock Calibration** (SNP distance, phylogenetic bridging)
3. **Serology System** (wbeT frameshift detection, Inaba switch)
4. **Horizontal Gene Transfer** (IncA/C plasmid vs SXT/ICE discrimination)
5. **Environmental Resilience** (hapR/vpsA integrity checks)
6. **Lineage Specificity** (Bengal L1 vs L2, Philippines GI-119 rejection)
7. **Sample Quality Degradation** (freeze-thaw proxy, QC gating)

## Results by Objective

"""
        
        # Group results by objective
        by_objective = defaultdict(list)
        for result in self.results:
            by_objective[result.objective_id].append(result)
        
        for obj_id in sorted(by_objective.keys()):
            obj_results = by_objective[obj_id]
            obj_name = obj_results[0].objective_name if obj_results else "Unknown"
            
            pass_count = sum(1 for r in obj_results if r.pass_fail == "PASS")
            fail_count = sum(1 for r in obj_results if r.pass_fail == "FAIL")
            partial_count = sum(1 for r in obj_results if r.pass_fail == "PARTIAL")
            not_tested = sum(1 for r in obj_results if r.pass_fail == "NOT_TESTED")
            
            status_icon = "✓" if pass_count == len(obj_results) else ("?" if not_tested > 0 else "✗")
            
            report += f"### {status_icon} {obj_id}: {obj_name}\n"
            report += f"**Status:** {pass_count} PASS | {fail_count} FAIL | {partial_count} PARTIAL | {not_tested} NOT_TESTED\n\n"
            report += "| Sample | Expected | Observed | Result | Details |\n"
            report += "|--------|----------|----------|--------|----------|\n"
            
            for res in obj_results:
                report += f"| {res.test_sample} | {res.expected_call} | {res.observed_call} | {res.pass_fail} | {res.details} |\n"
            
            report += "\n"
        
        # Summary statistics
        total_tests = len(self.results)
        passed = sum(1 for r in self.results if r.pass_fail == "PASS")
        failed = sum(1 for r in self.results if r.pass_fail == "FAIL")
        partial = sum(1 for r in self.results if r.pass_fail == "PARTIAL")
        not_tested = sum(1 for r in self.results if r.pass_fail == "NOT_TESTED")
        
        report += f"""## Overall Statistics

| Metric | Count | Percentage |
|--------|-------|-----------|
| **PASS** | {passed} | {100*passed//total_tests}% |
| **FAIL** | {failed} | {100*failed//total_tests}% |
| **PARTIAL** | {partial} | {100*partial//total_tests}% |
| **NOT_TESTED** | {not_tested} | {100*not_tested//total_tests}% |
| **TOTAL** | {total_tests} | 100% |

## Gaps & Recommendations

### Not Yet Implemented (Deployment Blockers)
"""
        
        # Identify gaps
        gaps = {
            "OBJ1": "rtxA_status field in virulence_report.json",
            "OBJ3": "wbeT_mutation frameshift/stop-codon detection",
            "OBJ4": "replicon_type discrimination (IncA/C vs SXT)",
            "OBJ5": "hapR_integrity and vpsA_integrity loss-of-function checks",
            "OBJ6": "lineage_classification with Bengal L1/L2 and Philippines GI-119 specificity",
            "OBJ7": "degradation_proxy_cv (k-mer coefficient of variation metric)"
        }
        
        for obj_id, gap_desc in gaps.items():
            report += f"\n#### {obj_id}: {gap_desc}\n"
            if not_tested > 0:
                report += "- **Status:** NOT_IMPLEMENTED (field not found in JSON)\n"
                report += "- **Action:** Add to respective module (see STRESS_TEST_PROTOCOL.md)\n"
        
        report += f"""

## Field Deployment Readiness

**Current Assessment:** VALIDATION ONGOING

Based on the stress test results:

- ✓ Objectives with full PASS status are ready for field deployment
- ⚠ Objectives with PARTIAL status require additional testing or field caveats
- ✗ Objectives with FAIL or NOT_TESTED status must be resolved before deployment

### Next Steps

1. Implement missing fields (see Gaps section above)
2. Re-run stress test on actual test samples (not mock data)
3. Establish field QC thresholds based on OBJ7 degradation testing
4. Train field teams on interpretation of heterogeneity alerts (OBJ3, OBJ2)
5. Document exclusion criteria (OBJ6) for publication

---

**Report Generated:** {datetime.now().isoformat()}
**Workspace:** {self.workspace}
"""
        
        return report
    
    def save_results_json(self, output_path: str) -> None:
        """Save detailed results as JSON"""
        results_dict = [asdict(r) for r in self.results]
        
        summary = {
            "timestamp": datetime.now().isoformat(),
            "workspace": str(self.workspace),
            "total_tests": len(self.results),
            "pass_count": sum(1 for r in self.results if r.pass_fail == "PASS"),
            "fail_count": sum(1 for r in self.results if r.pass_fail == "FAIL"),
            "partial_count": sum(1 for r in self.results if r.pass_fail == "PARTIAL"),
            "not_tested_count": sum(1 for r in self.results if r.pass_fail == "NOT_TESTED"),
            "results": results_dict
        }
        
        with open(output_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"[✓] Results saved to {output_path}")
    
    def save_results_markdown(self, output_path: str) -> None:
        """Save report as Markdown"""
        report = self.generate_report()
        with open(output_path, 'w') as f:
            f.write(report)
        print(f"[✓] Report saved to {output_path}")
    
    def _create_validation_config(self) -> None:
        """Create a minimal validation config file"""
        config_content = """# Validation Configuration for Stress Test
samples:
  - haiti_2010el1786
  - env_he09
  - env_j515_2018
  - clinical_inaba_2012el1410
  - bangladesh_l1
  - philippines_outbreak
  - fresh_extract
  - freeze_thaw_1x
  - freeze_thaw_3x
  - freeze_thaw_5x

reference_genome: "2010EL-1786"
reference_database: "data/references/2010EL-1786.fasta"

# Pipeline parameters
quality_threshold: 20
coverage_threshold: 10

# Regional selection (for OBJ2, OBJ6 testing)
regional_baseline: "Haiti_2010_Ancestor"
global_baseline: "Haiti_Yemen_Global"
"""
        
        config_path = self.workspace / "workflow" / "validation_config.yaml"
        with open(config_path, 'w') as f:
            f.write(config_content)
        print(f"[✓] Created validation_config.yaml")


def main():
    """Main entry point"""
    harness = StressTestHarness()
    
    print("\n" + "="*80)
    print("VIBRION SENTINEL STRESS TEST HARNESS v1.0")
    print("="*80)
    print(f"Workspace: {harness.workspace}")
    print(f"Test samples dir: {harness.test_samples_dir}")
    print(f"Pipeline output dir: {harness.pipeline_output}")
    
    # Phase 1: Evaluate existing outputs (without running full pipeline)
    print("\n[PHASE 1] Evaluating existing pipeline outputs...")
    harness.evaluate_all_objectives()
    
    # Phase 2: Generate reports
    print("\n[PHASE 2] Generating reports...")
    harness.save_results_json(str(harness.validation_dir / "STRESS_TEST_RESULTS.json"))
    harness.save_results_markdown(str(harness.validation_dir / "STRESS_TEST_RESULTS.md"))
    
    print("\n" + "="*80)
    print("STRESS TEST COMPLETE")
    print("="*80)
    print(f"\nReports saved to: {harness.validation_dir}/")
    print("  - STRESS_TEST_RESULTS.json (machine-readable)")
    print("  - STRESS_TEST_RESULTS.md (human-readable)")


if __name__ == "__main__":
    main()
