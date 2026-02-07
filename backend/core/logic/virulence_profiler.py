#!/usr/bin/env python3
"""
Enhanced virulence profiling for Vibrion cholerae
Detects specific toxin mutations including rtxA stop codons
"""

import json
from typing import Dict, Any, Optional
from pathlib import Path

class VirulenceProfiler:
    """Enhanced toxin and virulence marker detection"""
    
    # Key virulence mutations in V. cholerae
    RTXA_STOP_CODON_POSITION = 13602  # G13602A creates stop codon
    HEMOLYSIN_PRESENCE = "hlyA"
    TCP_PRESENCE = "tcpA"
    CTX_TOXIN = "ctxB"
    
    # Mutation signatures
    HAITI_RTXA_MUTATION = {
        "position": 13602,
        "ref_codon": "GAA",  # Glutamic acid
        "alt_codon": "TAA",  # STOP
        "effect": "Loss of function (MARTX toxin inactivation)",
        "lineage": "Haiti_2010_Ancestor"
    }
    
    def __init__(self):
        self.markers = {}
        
    def detect_rtxa_status(self, vcf_data: Dict[str, Any], reference_genome: Optional[str] = None) -> Dict[str, Any]:
        """
        Detect rtxA mutations, specifically the G13602A stop codon in Haiti strains
        
        Args:
            vcf_data: VCF variant calls or SNP dictionary
            reference_genome: FASTA file or sequence string
        
        Returns:
            Dictionary with rtxA status and mutation details
        """
        result = {
            "rtxA_status": "UNKNOWN",
            "rtxA_mutation": None,
            "rtxA_codon_position": self.RTXA_STOP_CODON_POSITION,
            "rtxA_effect": None,
            "is_haiti_ancestor": False,
            "confidence": "UNKNOWN"
        }
        
        # Check if we have variant calls at rtxA position
        if isinstance(vcf_data, dict) and "variants" in vcf_data:
            # Look for variant at position 13602
            for variant in vcf_data.get("variants", []):
                if variant.get("pos") == self.RTXA_STOP_CODON_POSITION:
                    ref = variant.get("ref", "")
                    alt = variant.get("alt", "")
                    
                    # Check for GAA→TAA mutation
                    if ref.upper() == "G" and alt.upper() == "A":
                        result["rtxA_status"] = "G13602A_STOP_CODON"
                        result["rtxA_mutation"] = "G13602A"
                        result["rtxA_codon_position"] = self.RTXA_STOP_CODON_POSITION
                        result["rtxA_effect"] = "MARTX toxin inactivation (premature stop codon)"
                        result["is_haiti_ancestor"] = True
                        result["confidence"] = "HIGH"
                        return result
        
        # Check for WT (wild-type) by absence of mutation
        # If we reach here without finding the mutation, assume WT
        # (This is a simplification; robust implementation would check coverage)
        result["rtxA_status"] = "FUNCTIONAL_WT"
        result["rtxA_mutation"] = "None (wild-type)"
        result["rtxA_effect"] = "MARTX toxin functional (standard El Tor)"
        result["is_haiti_ancestor"] = False
        result["confidence"] = "MEDIUM"  # Confidence lower for negative calls
        
        return result
    
    def detect_hemolysin_status(self, k_mer_matches: Dict[str, Any]) -> Dict[str, Any]:
        """Detect hemolysin (hlyA) presence and status"""
        result = {
            "hemolysin_status": "NOT_DETECTED",
            "hlyA_present": False,
            "hemolysin_function": None
        }
        
        if k_mer_matches.get(self.HEMOLYSIN_PRESENCE):
            result["hemolysin_status"] = "DETECTED"
            result["hlyA_present"] = True
            result["hemolysin_function"] = "Pore-forming toxin (alternative to MARTX)"
        
        return result
    
    def profile_virulence_strategy(self, rtxa_status: str, hemolysin_present: bool) -> Dict[str, Any]:
        """
        Classify virulence strategy based on toxin profile
        Haiti strains use hemolysin (due to MARTX inactivation)
        """
        result = {
            "virulence_strategy": "UNKNOWN",
            "threat_level": "UNKNOWN",
            "description": ""
        }
        
        if "STOP_CODON" in rtxa_status and hemolysin_present:
            result["virulence_strategy"] = "Haiti_Hypervirulent_Hemolysin_Dependent"
            result["threat_level"] = "HIGH"
            result["description"] = "MARTX-inactive strain relying on hemolysin; characteristic of Haiti 2010 epidemic strain"
        elif "STOP_CODON" in rtxa_status:
            result["virulence_strategy"] = "MARTX_Inactive_Alternative_Toxins_Unknown"
            result["threat_level"] = "MODERATE"
            result["description"] = "Possible Haiti-like variant; rtxA inactivated but hemolysin status uncertain"
        elif "FUNCTIONAL_WT" in rtxa_status and hemolysin_present:
            result["virulence_strategy"] = "Dual_Toxin_Standard_El_Tor"
            result["threat_level"] = "MODERATE"
            result["description"] = "Standard El Tor with both MARTX and hemolysin"
        elif "FUNCTIONAL_WT" in rtxa_status:
            result["virulence_strategy"] = "MARTX_Dependent_Standard_El_Tor"
            result["threat_level"] = "MODERATE"
            result["description"] = "Environmental or atypical strain with functional MARTX"
        
        return result


def test_rtxa_detection():
    """Unit test for rtxA detection logic"""
    profiler = VirulenceProfiler()
    
    # Test 1: Haiti strain (G13602A present)
    vcf_haiti = {
        "variants": [
            {"pos": 13602, "ref": "G", "alt": "A"},
        ]
    }
    result_haiti = profiler.detect_rtxa_status(vcf_haiti)
    assert result_haiti["rtxA_status"] == "G13602A_STOP_CODON"
    assert result_haiti["is_haiti_ancestor"] == True
    print("✓ Test 1 passed: Haiti G13602A detection")
    
    # Test 2: WT strain (no mutation)
    vcf_wt = {"variants": []}
    result_wt = profiler.detect_rtxa_status(vcf_wt)
    assert result_wt["rtxA_status"] == "FUNCTIONAL_WT"
    assert result_wt["is_haiti_ancestor"] == False
    print("✓ Test 2 passed: WT strain detection")
    
    # Test 3: Virulence strategy classification
    strategy = profiler.profile_virulence_strategy("G13602A_STOP_CODON", True)
    assert "Haiti_Hypervirulent" in strategy["virulence_strategy"]
    print("✓ Test 3 passed: Virulence strategy classification")


if __name__ == "__main__":
    test_rtxa_detection()
    print("\nAll tests passed!")
