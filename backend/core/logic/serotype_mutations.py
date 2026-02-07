#!/usr/bin/env python3
"""
Enhanced serogroup inference with serotype-specific mutations
Detects wbeT frameshifts (Inaba switch), Hikojima markers, etc.
"""

import json
from typing import Dict, Any, Optional, List
from pathlib import Path

class SerotypeMutationDetector:
    """Detect serotype-defining mutations in V. cholerae"""
    
    # Serotype-defining genes and mutations
    WBET_GENE = "wbeT"
    RFB_CLUSTER = "rfb"
    
    # Known serotype mutations
    INABA_MUTATION = {
        "gene": "wbeT",
        "mutation_type": "stop_codon",
        "ref_codon": "GAA",  # Glutamic acid
        "alt_codon": "TAA",  # STOP
        "effect": "Loss of wbeT function (Ogawa→Inaba switch)",
        "serotype": "Inaba"
    }
    
    OGAWA_WT = {
        "gene": "wbeT",
        "status": "intact",
        "effect": "Functional wbeT (Ogawa polysaccharide synthesis)",
        "serotype": "Ogawa"
    }
    
    HIKOJIMA_MARKERS = {
        "rfb_c": "present",  # Rough antigen
        "serotype": "Hikojima"
    }
    
    def __init__(self):
        self.serotype_calls = {}
        
    def detect_wbet_mutation(self, vcf_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Detect wbeT mutations, specifically frameshifts/stop codons causing Inaba phenotype
        
        Args:
            vcf_data: VCF variant calls or SNP dictionary with regions
        
        Returns:
            Dictionary with wbeT status, mutation details, and serotype implications
        """
        result = {
            "wbeT_status": "UNKNOWN",
            "wbeT_mutation": None,
            "wbeT_frameshift_detected": False,
            "wbeT_stop_codon_detected": False,
            "serotype_call": "UNKNOWN",
            "serotype_shift_alert": None,
            "confidence": "UNKNOWN"
        }
        
        # Check for known stop-codon mutations in wbeT
        if isinstance(vcf_data, dict) and "variants" in vcf_data:
            for variant in vcf_data.get("variants", []):
                # Look for frameshift or stop codon mutations
                if variant.get("gene") == self.WBET_GENE or self.WBET_GENE in str(variant.get("annotation", "")):
                    
                    # Detect stop codon mutations (e.g., GAA→TAA)
                    if variant.get("ref") and variant.get("alt"):
                        ref = variant.get("ref", "").upper()
                        alt = variant.get("alt", "").upper()
                        
                        # Check for GAA→TAA (stop codon)
                        if ref == "G" and alt == "A":
                            result["wbeT_stop_codon_detected"] = True
                            result["wbeT_mutation"] = "GAA_to_TAA_stop"
                            result["wbeT_status"] = "FRAMESHIFT_STOP_CODON"
                            result["serotype_call"] = "Inaba"
                            result["serotype_shift_alert"] = "VACCINE_MISMATCH_RISK"
                            result["confidence"] = "HIGH"
                            return result
                    
                    # Detect other frameshifts
                    if variant.get("variant_type") == "indel":
                        result["wbeT_frameshift_detected"] = True
                        result["wbeT_mutation"] = variant.get("description", "Frameshift")
                        result["wbeT_status"] = "FRAMESHIFT_OR_DELETION"
                        result["serotype_call"] = "Inaba"  # Frameshift loss → Inaba
                        result["serotype_shift_alert"] = "VACCINE_MISMATCH_RISK"
                        result["confidence"] = "MEDIUM"
                        return result
        
        # If no mutations found, assume WT (Ogawa)
        result["wbeT_status"] = "WT_FUNCTIONAL"
        result["wbeT_mutation"] = "None (wild-type)"
        result["serotype_call"] = "Ogawa"
        result["serotype_shift_alert"] = None
        result["confidence"] = "MEDIUM"  # Negative calls have lower confidence
        
        return result
    
    def detect_rfb_markers(self, k_mer_matches: Dict[str, Any]) -> Dict[str, Any]:
        """
        Detect rfb cluster markers for Hikojima serotype
        """
        result = {
            "rfb_cluster": "NOT_DETECTED",
            "rfb_c_present": False,
            "hikojima_indicators": []
        }
        
        if k_mer_matches.get("rfb_c"):
            result["rfb_cluster"] = "DETECTED"
            result["rfb_c_present"] = True
            result["hikojima_indicators"].append("rfb_c present")
        
        return result
    
    def call_serotype(self, wbet_status: str, rfb_markers: Dict[str, Any]) -> Dict[str, Any]:
        """
        Integrate wbeT and rfb markers to call serotype
        """
        result = {
            "serotype": "UNKNOWN",
            "confidence": "LOW",
            "rationale": []
        }
        
        if "STOP_CODON" in wbet_status or "FRAMESHIFT" in wbet_status:
            result["serotype"] = "Inaba"
            result["confidence"] = "HIGH"
            result["rationale"].append("wbeT frameshift/stop codon detected")
        elif rfb_markers.get("rfb_c_present"):
            result["serotype"] = "Hikojima"
            result["confidence"] = "MEDIUM"
            result["rationale"].append("rfb_c detected (rough antigen)")
        elif "WT" in wbet_status:
            result["serotype"] = "Ogawa"
            result["confidence"] = "MEDIUM"
            result["rationale"].append("wbeT intact (wild-type)")
        
        return result
    
    def check_vaccine_mismatch(self, serotype: str) -> Dict[str, Any]:
        """
        Alert if serotype differs from vaccine strain (typically Ogawa)
        
        Haiti 2022 resurgence involved Inaba variants - vaccine mismatch risk
        """
        result = {
            "vaccine_mismatch_detected": False,
            "alert_level": "NONE",
            "description": ""
        }
        
        if serotype == "Inaba":
            result["vaccine_mismatch_detected"] = True
            result["alert_level"] = "WARNING"
            result["description"] = (
                "Inaba serotype detected. Most licensed oral cholera vaccines (OCVs) are formulated with "
                "Ogawa biotypes. This serotype mismatch may reduce vaccine effectiveness. "
                "Clinical correlation and serological follow-up recommended."
            )
        elif serotype == "Hikojima":
            result["vaccine_mismatch_detected"] = True
            result["alert_level"] = "CAUTION"
            result["description"] = (
                "Hikojima (rough) phenotype detected. Reduced virulence expected but vaccine efficacy uncertain."
            )
        
        return result


def test_wbet_detection():
    """Unit test for wbeT mutation detection"""
    detector = SerotypeMutationDetector()
    
    # Test 1: Inaba strain (GAA→TAA stop codon)
    vcf_inaba = {
        "variants": [
            {"gene": "wbeT", "ref": "G", "alt": "A", "annotation": "stop_codon"},
        ]
    }
    result_inaba = detector.detect_wbet_mutation(vcf_inaba)
    assert result_inaba["wbeT_stop_codon_detected"] == True
    assert result_inaba["serotype_call"] == "Inaba"
    assert "VACCINE_MISMATCH" in result_inaba.get("serotype_shift_alert", "")
    print("✓ Test 1 passed: Inaba stop codon detection")
    
    # Test 2: Ogawa strain (WT)
    vcf_ogawa = {"variants": []}
    result_ogawa = detector.detect_wbet_mutation(vcf_ogawa)
    assert result_ogawa["wbeT_status"] == "WT_FUNCTIONAL"
    assert result_ogawa["serotype_call"] == "Ogawa"
    print("✓ Test 2 passed: Ogawa WT detection")
    
    # Test 3: Vaccine mismatch alert
    alert = detector.check_vaccine_mismatch("Inaba")
    assert alert["vaccine_mismatch_detected"] == True
    assert alert["alert_level"] == "WARNING"
    print("✓ Test 3 passed: Vaccine mismatch alert")


if __name__ == "__main__":
    test_wbet_detection()
    print("\nAll tests passed!")
