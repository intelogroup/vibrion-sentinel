#!/usr/bin/env python3
"""
Integration module: Connect all 7 stress test validation modules to the comprehensive reporter
Extends generate_comprehensive_report.py with new biological checks
"""

import json
import sys
from pathlib import Path
from typing import Dict, Any, Optional

# Import all validation modules (relative imports for package context)
try:
    # Try package imports first (when run as module)
    from backend.core.logic.virulence_profiler import VirulenceProfiler
    from backend.core.logic.serotype_mutations import SerotypeMutationDetector
    from backend.core.logic.amr_element_discriminator import AMRElementDiscriminator
    from backend.core.logic.environmental_resilience import EnvironmentalResilienceProfiler
    from backend.core.logic.lineage_specificity import LineageSpecificityClassifier
    from backend.core.logic.degradation_proxy import DegradationProxyCalculator, SNPDistanceCalculator
except ImportError:
    try:
        # Try direct imports (when run as script)
        from virulence_profiler import VirulenceProfiler
        from serotype_mutations import SerotypeMutationDetector
        from amr_element_discriminator import AMRElementDiscriminator
        from environmental_resilience import EnvironmentalResilienceProfiler
        from lineage_specificity import LineageSpecificityClassifier
        from degradation_proxy import DegradationProxyCalculator, SNPDistanceCalculator
    except ImportError as e:
        print(f"Warning: Could not import validation modules: {e}")
        print("Stress test functionality will be limited")
        VirulenceProfiler = None
        SerotypeMutationDetector = None
        AMRElementDiscriminator = None
        EnvironmentalResilienceProfiler = None
        LineageSpecificityClassifier = None
        DegradationProxyCalculator = None
        SNPDistanceCalculator = None


class StressTestIntegrator:
    """Integrate all 7 validation objectives into report generation"""
    
    def __init__(self):
        self.profilers = {}
        
        # Initialize profilers only if modules loaded successfully
        if VirulenceProfiler:
            self.profilers["virulence"] = VirulenceProfiler()
        if SerotypeMutationDetector:
            self.profilers["serotype"] = SerotypeMutationDetector()
        if AMRElementDiscriminator:
            self.profilers["amr"] = AMRElementDiscriminator()
        if EnvironmentalResilienceProfiler:
            self.profilers["resilience"] = EnvironmentalResilienceProfiler()
        if LineageSpecificityClassifier:
            self.profilers["lineage"] = LineageSpecificityClassifier()
        if DegradationProxyCalculator:
            self.profilers["qc"] = DegradationProxyCalculator()
        if SNPDistanceCalculator:
            self.profilers["snp_distance"] = SNPDistanceCalculator()
    
    def run_all_validations(self, vcf_data: Dict[str, Any], k_mer_matches: Dict[str, Any],
                           vibrio_stats: Dict[str, Any], kmer_depths: Optional[list] = None) -> Dict[str, Any]:
        """
        Execute all 7 biological validation objectives
        """
        results = {
            "timestamp": None,
            "objective_1_rtxa": {},
            "objective_2_clock": {},
            "objective_3_serotype": {},
            "objective_4_hgt": {},
            "objective_5_resilience": {},
            "objective_6_lineage": {},
            "objective_7_qc": {}
        }
        
        # OBJ1: rtxA Stop-Codon Detection
        try:
            rtxa_result = self.profilers["virulence"].detect_rtxa_status(vcf_data)
            hemolysin = self.profilers["virulence"].detect_hemolysin_status(k_mer_matches)
            virulence_strategy = self.profilers["virulence"].profile_virulence_strategy(
                rtxa_result["rtxA_status"], hemolysin["hlyA_present"]
            )
            results["objective_1_rtxa"] = {
                "rtxa_status": rtxa_result,
                "hemolysin_status": hemolysin,
                "virulence_strategy": virulence_strategy
            }
        except Exception as e:
            results["objective_1_rtxa"]["error"] = str(e)
        
        # OBJ2: SNP Distance & Time Calibration
        try:
            snp_dist = self.profilers["snp_distance"].calculate_snp_distance(vcf_data)
            divergence = self.profilers["snp_distance"].estimate_divergence_date(snp_dist)
            results["objective_2_clock"] = {
                "snp_distance_to_reference": snp_dist,
                "divergence_estimate": divergence,
                "phylo_position": self._infer_phylo_position(snp_dist)
            }
        except Exception as e:
            results["objective_2_clock"]["error"] = str(e)
        
        # OBJ3: wbeT Frameshift & Serotype
        try:
            wbet_result = self.profilers["serotype"].detect_wbet_mutation(vcf_data)
            rfb_markers = self.profilers["serotype"].detect_rfb_markers(k_mer_matches)
            serotype_call = self.profilers["serotype"].call_serotype(wbet_result["wbeT_status"], rfb_markers)
            vaccine_mismatch = self.profilers["serotype"].check_vaccine_mismatch(serotype_call["serotype"])
            results["objective_3_serotype"] = {
                "wbet_status": wbet_result,
                "rfb_markers": rfb_markers,
                "serotype_call": serotype_call,
                "vaccine_mismatch": vaccine_mismatch
            }
        except Exception as e:
            results["objective_3_serotype"]["error"] = str(e)
        
        # OBJ4: IncA/C vs SXT/ICE Discrimination
        try:
            sxt_result = self.profilers["amr"].detect_sxt_element({}, k_mer_matches)
            inca_result = self.profilers["amr"].detect_inca_c_plasmid(k_mer_matches)
            amr_class = self.profilers["amr"].classify_amr_element(sxt_result, inca_result)
            transmission = self.profilers["amr"].predict_transmission_dynamics(amr_class["element_type"])
            results["objective_4_hgt"] = {
                "sxt_element": sxt_result,
                "inca_c_plasmid": inca_result,
                "classification": amr_class,
                "transmission_dynamics": transmission
            }
        except Exception as e:
            results["objective_4_hgt"]["error"] = str(e)
        
        # OBJ5: hapR/vpsA Integrity & Biofilm
        try:
            hapr_result = self.profilers["resilience"].check_hapr_integrity(vcf_data)
            vpsa_result = self.profilers["resilience"].check_vpsa_integrity(vcf_data)
            vps_cluster = self.profilers["resilience"].check_vps_cluster_completeness(k_mer_matches)
            resilience_report = self.profilers["resilience"].generate_resilience_report(
                hapr_result, vpsa_result, vps_cluster
            )
            results["objective_5_resilience"] = {
                "hapR_integrity": hapr_result,
                "vpsA_integrity": vpsa_result,
                "vps_cluster": vps_cluster,
                "biofilm_profile": resilience_report
            }
        except Exception as e:
            results["objective_5_resilience"]["error"] = str(e)
        
        # OBJ6: Lineage Specificity & Imposter Detection
        try:
            snp_dist_for_lineage = results.get("objective_2_clock", {}).get("snp_distance_to_reference")
            lineage_result = self.profilers["lineage"].classify_lineage(vcf_data, k_mer_matches, snp_dist_for_lineage)
            lineage_report = self.profilers["lineage"].generate_lineage_report(lineage_result)
            results["objective_6_lineage"] = {
                "classification": lineage_result,
                "report": lineage_report
            }
        except Exception as e:
            results["objective_6_lineage"]["error"] = str(e)
        
        # OBJ7: Sample Quality & Degradation Proxy
        try:
            qc_report = self.profilers["qc"].generate_qc_report(vibrio_stats, {}, kmer_depths)
            results["objective_7_qc"] = qc_report
        except Exception as e:
            results["objective_7_qc"]["error"] = str(e)
        
        return results
    
    def _infer_phylo_position(self, snp_distance: float) -> str:
        """Infer phylogenetic position from SNP distance"""
        if snp_distance == 0:
            return "ROOT_ANCESTOR (2010 Haiti)"
        elif snp_distance < 10:
            return "VERY_BASAL (early divergence)"
        elif snp_distance < 20:
            return "BASAL_NODE_BRIDGE (2015-2018 environmental bridge)"
        elif snp_distance < 35:
            return "INTERMEDIATE (2018-2020 transition)"
        elif snp_distance < 60:
            return "DERIVED (late epidemic strains)"
        else:
            return "HIGHLY_DIVERGED (possible foreign lineage)"


def test_integrator():
    """Unit test for stress test integrator"""
    integrator = StressTestIntegrator()
    
    # Mock data
    vcf_data = {
        "variants": [
            {"pos": 13602, "ref": "G", "alt": "A"}  # rtxA mutation
        ]
    }
    
    k_mer_matches = {
        "setCD": True,
        "acaCD": True,
        "floR": True,
        "hlyA": True
    }
    
    vibrio_stats = {
        "vibrio_percentage": 98.5,
        "total_reads": 1000000
    }
    
    results = integrator.run_all_validations(vcf_data, k_mer_matches, vibrio_stats)
    
    print(json.dumps(results, indent=2, default=str))
    
    # Verify results
    assert "objective_1_rtxa" in results
    assert results["objective_1_rtxa"]["rtxa_status"]["rtxA_status"] == "G13602A_STOP_CODON"
    print("✓ Integration test passed")


if __name__ == "__main__":
    test_integrator()
