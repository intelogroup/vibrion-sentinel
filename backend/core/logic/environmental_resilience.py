#!/usr/bin/env python3
"""
Environmental resilience profiling: hapR and vpsA integrity checks
Detects loss-of-function mutations affecting biofilm phenotype
"""

import json
from typing import Dict, Any, Optional, List
from pathlib import Path

class EnvironmentalResilienceProfiler:
    """Detect functional integrity of biofilm/environmental adaptation genes"""
    
    # Key environmental genes
    HAPR_GENE = "hapR"  # Quorum sensing regulator
    VPSA_GENE = "vpsA"  # Biofilm polysaccharide biosynthesis
    VPS_CLUSTER = ["vpsA", "vpsB", "vpsC", "vpsD", "vpsE"]
    
    def __init__(self):
        self.biofilm_profile = {}
        
    def check_hapr_integrity(self, vcf_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Check for loss-of-function mutations in hapR (quorum sensing regulator)
        HapR is essential for biofilm regulation in V. cholerae
        """
        result = {
            "hapR_integrity": "UNKNOWN",
            "hapR_status": "UNKNOWN",
            "hapR_mutation_type": None,
            "hapR_has_frameshift": False,
            "hapR_has_stop_codon": False,
            "hapR_functional": None,
            "confidence": "UNKNOWN"
        }
        
        if isinstance(vcf_data, dict) and "variants" in vcf_data:
            lof_mutations = 0
            
            for variant in vcf_data.get("variants", []):
                if variant.get("gene") == self.HAPR_GENE or self.HAPR_GENE in str(variant.get("annotation", "")):
                    
                    # Check for stop codons (nonsense mutations)
                    if variant.get("variant_type") == "stop_codon" or "stop" in str(variant.get("effect", "")).lower():
                        result["hapR_has_stop_codon"] = True
                        result["hapR_mutation_type"] = "Stop codon (nonsense)"
                        lof_mutations += 1
                    
                    # Check for frameshifts
                    if variant.get("variant_type") == "frameshift" or variant.get("variant_type") == "indel":
                        if variant.get("size", 0) % 3 != 0:  # Not multiple of 3
                            result["hapR_has_frameshift"] = True
                            result["hapR_mutation_type"] = "Frameshift"
                            lof_mutations += 1
            
            if lof_mutations > 0:
                result["hapR_integrity"] = "LOF_STOP_CODON" if result["hapR_has_stop_codon"] else "LOF_FRAMESHIFT"
                result["hapR_status"] = "Loss of function"
                result["hapR_functional"] = False
                result["confidence"] = "HIGH"
            else:
                result["hapR_integrity"] = "WT_PROFICIENT"
                result["hapR_status"] = "Wild-type (functional)"
                result["hapR_functional"] = True
                result["confidence"] = "MEDIUM"
        else:
            # No variant data; assume WT
            result["hapR_integrity"] = "WT_PROFICIENT"
            result["hapR_status"] = "Assumed wild-type (no mutations detected)"
            result["hapR_functional"] = True
            result["confidence"] = "LOW"
        
        return result
    
    def check_vpsa_integrity(self, vcf_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Check for loss-of-function mutations in vpsA (biofilm matrix synthesis)
        VpsA is essential for the rugose phenotype
        """
        result = {
            "vpsA_integrity": "UNKNOWN",
            "vpsA_status": "UNKNOWN",
            "vpsA_mutation_type": None,
            "vpsA_has_frameshift": False,
            "vpsA_has_stop_codon": False,
            "vpsA_functional": None,
            "confidence": "UNKNOWN"
        }
        
        if isinstance(vcf_data, dict) and "variants" in vcf_data:
            lof_mutations = 0
            
            for variant in vcf_data.get("variants", []):
                if variant.get("gene") == self.VPSA_GENE or self.VPSA_GENE in str(variant.get("annotation", "")):
                    
                    # Check for stop codons
                    if variant.get("variant_type") == "stop_codon" or "stop" in str(variant.get("effect", "")).lower():
                        result["vpsA_has_stop_codon"] = True
                        result["vpsA_mutation_type"] = "Stop codon (nonsense)"
                        lof_mutations += 1
                    
                    # Check for frameshifts
                    if variant.get("variant_type") == "frameshift" or variant.get("variant_type") == "indel":
                        if variant.get("size", 0) % 3 != 0:
                            result["vpsA_has_frameshift"] = True
                            result["vpsA_mutation_type"] = "Frameshift"
                            lof_mutations += 1
            
            if lof_mutations > 0:
                result["vpsA_integrity"] = "LOF_STOP_CODON" if result["vpsA_has_stop_codon"] else "LOF_FRAMESHIFT"
                result["vpsA_status"] = "Loss of function"
                result["vpsA_functional"] = False
                result["confidence"] = "HIGH"
            else:
                result["vpsA_integrity"] = "WT_FUNCTIONAL"
                result["vpsA_status"] = "Wild-type (functional)"
                result["vpsA_functional"] = True
                result["confidence"] = "MEDIUM"
        else:
            result["vpsA_integrity"] = "WT_FUNCTIONAL"
            result["vpsA_status"] = "Assumed wild-type"
            result["vpsA_functional"] = True
            result["confidence"] = "LOW"
        
        return result
    
    def check_vps_cluster_completeness(self, k_mer_matches: Dict[str, Any]) -> Dict[str, Any]:
        """
        Check presence of the entire vps cluster for biofilm synthesis
        """
        result = {
            "vps_cluster_detected": False,
            "vps_genes_found": [],
            "vps_genes_missing": [],
            "cluster_completeness": 0.0
        }
        
        if isinstance(k_mer_matches, dict):
            for gene in self.VPS_CLUSTER:
                if k_mer_matches.get(gene):
                    result["vps_genes_found"].append(gene)
                else:
                    result["vps_genes_missing"].append(gene)
            
            completeness = len(result["vps_genes_found"]) / len(self.VPS_CLUSTER)
            result["cluster_completeness"] = completeness
            result["vps_cluster_detected"] = completeness > 0.7  # >70% detected
        
        return result
    
    def predict_biofilm_phenotype(self, hapr_functional: bool, vpsa_functional: bool, vps_complete: bool) -> Dict[str, Any]:
        """
        Predict biofilm phenotype based on regulatory and structural integrity
        
        Haiti 2010 strain: HapR proficient, VpsA intact → Rugose phenotype
        """
        result = {
            "biofilm_phenotype": "UNKNOWN",
            "environmental_resilience": "UNKNOWN",
            "growth_characteristics": [],
            "survival_prediction": "UNKNOWN"
        }
        
        if hapr_functional and vpsa_functional and vps_complete:
            result["biofilm_phenotype"] = "Rugose"
            result["environmental_resilience"] = "HIGH"
            result["growth_characteristics"] = [
                "Temperature-dependent: favored at 37°C",
                "Biofilm formation in saltwater",
                "Enhanced survival in environmental reservoirs"
            ]
            result["survival_prediction"] = "Extended survival in aquatic environment"
        elif not hapr_functional or not vpsa_functional:
            result["biofilm_phenotype"] = "Smooth"
            result["environmental_resilience"] = "REDUCED"
            result["growth_characteristics"] = [
                "Loss of quorum sensing regulation (hapR LoF)",
                "Reduced capsule/matrix synthesis (vpsA LoF)",
                "Planktonic phenotype predicted"
            ]
            result["survival_prediction"] = "Short-term survival; poor environmental persistence"
        elif not vps_complete:
            result["biofilm_phenotype"] = "Intermediate"
            result["environmental_resilience"] = "MODERATE"
            result["growth_characteristics"] = ["Incomplete vps cluster"]
            result["survival_prediction"] = "Reduced biofilm capacity"
        
        return result
    
    def generate_resilience_report(self, hapr_result: Dict, vpsa_result: Dict, vps_result: Dict) -> Dict[str, Any]:
        """
        Comprehensive environmental resilience assessment
        """
        hapr_func = hapr_result.get("hapR_functional", False)
        vpsa_func = vpsa_result.get("vpsA_functional", False)
        vps_complete = vps_result.get("cluster_completeness", 0.0) > 0.7
        
        biofilm_pred = self.predict_biofilm_phenotype(hapr_func, vpsa_func, vps_complete)
        
        report = {
            "hapR_integrity": hapr_result.get("hapR_integrity", "UNKNOWN"),
            "vpsA_integrity": vpsa_result.get("vpsA_integrity", "UNKNOWN"),
            "vps_cluster_completeness": vps_result.get("cluster_completeness", 0.0),
            "predicted_biofilm_phenotype": biofilm_pred.get("biofilm_phenotype"),
            "environmental_resilience": biofilm_pred.get("environmental_resilience"),
            "survival_prediction": biofilm_pred.get("survival_prediction"),
            "functional_assessment": {
                "hapR_proficient": hapr_func,
                "vpsA_functional": vpsa_func,
                "vps_pathway_intact": vps_complete
            }
        }
        
        return report


def test_resilience_detection():
    """Unit test for environmental resilience checks"""
    profiler = EnvironmentalResilienceProfiler()
    
    # Test 1: WT strain (Haiti-like)
    vcf_wt = {"variants": []}
    hapr_wt = profiler.check_hapr_integrity(vcf_wt)
    vpsa_wt = profiler.check_vpsa_integrity(vcf_wt)
    
    assert hapr_wt["hapR_functional"] == True
    assert vpsa_wt["vpsA_functional"] == True
    print("✓ Test 1 passed: WT strain (Haiti-like) detection")
    
    # Test 2: hapR LoF
    vcf_hapR_lof = {
        "variants": [
            {"gene": "hapR", "variant_type": "stop_codon", "effect": "stop_gained"}
        ]
    }
    hapr_lof = profiler.check_hapr_integrity(vcf_hapR_lof)
    assert hapr_lof["hapR_functional"] == False
    assert "LOF" in hapr_lof["hapR_integrity"]
    print("✓ Test 2 passed: hapR LoF detection")
    
    # Test 3: Biofilm phenotype prediction
    phenotype = profiler.predict_biofilm_phenotype(True, True, True)
    assert phenotype["biofilm_phenotype"] == "Rugose"
    assert phenotype["environmental_resilience"] == "HIGH"
    print("✓ Test 3 passed: Rugose phenotype prediction")
    
    # Test 4: Smooth phenotype (hapR LoF)
    phenotype_smooth = profiler.predict_biofilm_phenotype(False, True, True)
    assert phenotype_smooth["biofilm_phenotype"] == "Smooth"
    assert phenotype_smooth["environmental_resilience"] == "REDUCED"
    print("✓ Test 4 passed: Smooth phenotype prediction")


if __name__ == "__main__":
    test_resilience_detection()
    print("\nAll tests passed!")
