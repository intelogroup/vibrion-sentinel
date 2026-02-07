#!/usr/bin/env python3
"""
Lineage specificity: Distinguish Haiti/Global lineage from foreign variants
Includes Bengal L1 vs L2 and Philippines GI-119 discriminators
"""

import json
from typing import Dict, Any, Optional, List
from pathlib import Path

class LineageSpecificityClassifier:
    """Classify V. cholerae lineages and reject foreign imposters"""
    
    # Lineage-specific markers
    HAITI_GLOBAL_MARKERS = {
        "lineage_name": "Haiti_Yemen_Global",
        "lineage_id": "L2",  # Lineage 2
        "distinctive_features": [
            "Altered El Tor (ctxB7, rstR mutant)",
            "SXT/R391 ICE present",
            "VNTR profile: Haiti-specific haplotype",
            "SNP distance to 2010 ancestor: 0-50 SNPs"
        ],
        "outbreak_region": "Haiti 2010-present, Yemen 2016+, Global pandemic"
    }
    
    BENGAL_L1_MARKERS = {
        "lineage_name": "Bengal_L1_Dhaka_Endemic",
        "lineage_id": "L1",
        "distinctive_features": [
            "Tetracycline resistance via tetA(D) (not floR)",
            "Unique restriction-modification system (not GI-119)",
            "Distinct VNTR profile from Haiti",
            "SNP distance to Haiti: 50-100+ SNPs",
            "Geographic clustering in Dhaka, Bangladesh"
        ],
        "outbreak_region": "Bangladesh (endemic Dhaka population)"
    }
    
    PHILIPPINES_MARKERS = {
        "lineage_name": "Philippines_Hybrid_El_Tor",
        "lineage_id": "PHIL",
        "distinctive_features": [
            "Hybrid El Tor phenotype",
            "LACKS SXT element (unique)",
            "Carries genomic island GI-119 (diagnostic)",
            "Unique plasmid profile (IncA/C or other)",
            "SNP distance to Haiti: >100 SNPs",
            "Geographic clustering in Philippines"
        ],
        "outbreak_region": "Philippines (isolated outbreak)"
    }
    
    def __init__(self):
        self.lineage_calls = {}
        
    def detect_gi119(self, k_mer_matches: Dict[str, Any]) -> bool:
        """
        Detect GI-119 (Genomic Island 119) specific to Philippines strains
        """
        gi119_markers = ["GI119", "gi119", "genomic_island_119", "phi_119"]
        
        if isinstance(k_mer_matches, dict):
            for marker in gi119_markers:
                if k_mer_matches.get(marker):
                    return True
        
        return False
    
    def check_sxt_presence(self, k_mer_matches: Dict[str, Any]) -> bool:
        """Check if SXT element is present (should be absent in Philippines)"""
        sxt_markers = ["sxt", "setCD", "acaCD", "xis", "int"]
        
        if isinstance(k_mer_matches, dict):
            for marker in sxt_markers:
                if k_mer_matches.get(marker):
                    return True
        
        return False
    
    def detect_bengal_l1_specific_markers(self, k_mer_matches: Dict[str, Any]) -> Dict[str, bool]:
        """
        Detect markers specific to Bengal L1 (vs L2/Haiti)
        L1 has distinct R-M systems and tetA(D) preference
        """
        result = {
            "bengal_l1_rm_system": False,
            "tetA_D_only": False,
            "flo_R_absent": False
        }
        
        if isinstance(k_mer_matches, dict):
            # L1-specific restriction-modification system
            if k_mer_matches.get("bengal_l1_rm_system"):
                result["bengal_l1_rm_system"] = True
            
            # L1 uses tetA(D), not floR
            if k_mer_matches.get("tetA(D)") and not k_mer_matches.get("floR"):
                result["tetA_D_only"] = True
                result["flo_R_absent"] = True
        
        return result
    
    def classify_lineage(self, vcf_data: Dict[str, Any], k_mer_matches: Dict[str, Any], 
                        snp_distance: Optional[float] = None, phylo_context: Optional[str] = None) -> Dict[str, Any]:
        """
        Classify lineage based on markers, SNP distance, and phylogenetic context
        """
        result = {
            "lineage_classification": "UNKNOWN",
            "lineage_id": "UNKNOWN",
            "acceptance_status": "UNKNOWN",
            "confidence": "LOW",
            "supporting_evidence": [],
            "rejection_reason": None,
            "gi119_detected": False,
            "sxt_present": False
        }
        
        # Check for GI-119 (Philippines)
        result["gi119_detected"] = self.detect_gi119(k_mer_matches)
        result["sxt_present"] = self.check_sxt_presence(k_mer_matches)
        
        # Rule 1: If GI-119 present → Philippines
        if result["gi119_detected"]:
            result["lineage_classification"] = "Philippines_Hybrid_El_Tor"
            result["lineage_id"] = "PHIL"
            result["acceptance_status"] = "REJECT_FOREIGN"
            result["rejection_reason"] = "GI-119 marker detected; non-endemic Philippines lineage"
            result["confidence"] = "HIGH"
            result["supporting_evidence"].append("GI-119 genomic island signature")
            return result
        
        # Rule 2: If SXT absent AND GI-119 absent → possibly Philippines (but less certain)
        if not result["sxt_present"] and not result["gi119_detected"]:
            result["rejection_reason"] = "SXT element absent (unusual for Haiti/Bangladesh strains)"
            result["confidence"] = "MEDIUM"
        
        # Rule 3: Check for Bengal L1 markers
        bengal_markers = self.detect_bengal_l1_specific_markers(k_mer_matches)
        if bengal_markers["bengal_l1_rm_system"] or bengal_markers["tetA_D_only"]:
            result["lineage_classification"] = "Bengal_L1_Dhaka_Endemic"
            result["lineage_id"] = "L1"
            result["acceptance_status"] = "REJECT_FOREIGN"
            result["rejection_reason"] = "Bengal L1-specific markers detected; non-global-pandemic lineage"
            result["confidence"] = "MEDIUM_HIGH"
            result["supporting_evidence"].append("Bengal L1-specific R-M system")
            result["supporting_evidence"].append("tetA(D) preference (L1 profile)")
            return result
        
        # Rule 4: Default to Haiti/Global L2
        if snp_distance is not None and snp_distance < 100:
            result["lineage_classification"] = "Haiti_Yemen_Global_L2"
            result["lineage_id"] = "L2"
            result["acceptance_status"] = "ACCEPT_ENDEMIC"
            result["confidence"] = "HIGH"
            result["supporting_evidence"].append(f"SNP distance to 2010: {snp_distance:.1f} (< 100)")
            result["supporting_evidence"].append("SXT element present (Haiti/Global signature)")
        else:
            result["lineage_classification"] = "Haiti_Yemen_Global_L2"
            result["lineage_id"] = "L2"
            result["acceptance_status"] = "ACCEPT_ENDEMIC"
            result["confidence"] = "MEDIUM"
            result["supporting_evidence"].append("SXT element present; consistent with Haiti lineage")
        
        return result
    
    def generate_lineage_report(self, classification: Dict[str, Any]) -> Dict[str, Any]:
        """
        Generate detailed lineage report with clinical implications
        """
        lineage = classification.get("lineage_classification", "UNKNOWN")
        
        report = {
            "lineage": lineage,
            "status": classification.get("acceptance_status"),
            "clinical_interpretation": {},
            "field_action": None
        }
        
        if "L2" in lineage or "Haiti" in lineage or "Global" in lineage:
            report["clinical_interpretation"] = {
                "outbreak_context": "2010 Haiti epidemic lineage, currently endemic in Haiti/Yemen and circulating globally",
                "virulence_profile": "High (Classic toxin profile)",
                "vaccine_applicability": "Ogawa vaccines may not match if serotype is Inaba",
                "antimicrobial_profile": "SXT-mediated resistance typical"
            }
            report["field_action"] = "ACCEPT - Endemic strain, standard surveillance protocols apply"
        
        elif "L1" in lineage or "Bengal" in lineage:
            report["clinical_interpretation"] = {
                "outbreak_context": "Bangladesh-endemic lineage (non-global pandemic)",
                "virulence_profile": "Standard El Tor",
                "vaccine_applicability": "Different serotype profile from global pandemic",
                "antimicrobial_profile": "tetA(D) tetracycline resistance; distinct from SXT"
            }
            report["field_action"] = "REJECT - Non-endemic foreign lineage; distinct epidemiology from Haiti endemic"
        
        elif "Philippines" in lineage:
            report["clinical_interpretation"] = {
                "outbreak_context": "Geographically isolated Philippines outbreak",
                "virulence_profile": "Hybrid El Tor",
                "vaccine_applicability": "Unknown; distinct from global pandemic lineage",
                "antimicrobial_profile": "Lacks SXT; unique plasmid profile"
            }
            report["field_action"] = "REJECT - Foreign lineage; GI-119 signature confirms Philippines origin"
        
        else:
            report["field_action"] = "MANUAL_REVIEW - Lineage uncertain"
        
        return report


def test_lineage_classification():
    """Unit test for lineage specificity"""
    classifier = LineageSpecificityClassifier()
    
    # Test 1: Haiti/Global L2 classification
    vcf_haiti = {}
    k_mer_haiti = {"setCD": True, "acaCD": True, "floR": True}
    result_haiti = classifier.classify_lineage(vcf_haiti, k_mer_haiti, snp_distance=25)
    
    assert result_haiti["lineage_id"] == "L2"
    assert result_haiti["acceptance_status"] == "ACCEPT_ENDEMIC"
    print("✓ Test 1 passed: Haiti L2 classification")
    
    # Test 2: Philippines with GI-119
    k_mer_phil = {"GI119": True}
    result_phil = classifier.classify_lineage({}, k_mer_phil)
    
    assert result_phil["gi119_detected"] == True
    assert result_phil["lineage_id"] == "PHIL"
    assert result_phil["acceptance_status"] == "REJECT_FOREIGN"
    print("✓ Test 2 passed: Philippines rejection (GI-119)")
    
    # Test 3: Bengal L1 classification
    k_mer_bengal_l1 = {"bengal_l1_rm_system": True, "tetA(D)": True}
    result_l1 = classifier.classify_lineage({}, k_mer_bengal_l1)
    
    assert result_l1["lineage_id"] == "L1"
    assert result_l1["acceptance_status"] == "REJECT_FOREIGN"
    print("✓ Test 3 passed: Bengal L1 rejection")


if __name__ == "__main__":
    test_lineage_classification()
    print("\nAll tests passed!")
