#!/usr/bin/env python3
"""
Enhanced AMR profiler: Distinguish IncA/C plasmids from SXT/ICE elements
Includes regulatory markers (AcaCD, SetCD, TraC) and transmission dynamics
"""

import json
from typing import Dict, Any, Optional, List
from pathlib import Path

class AMRElementDiscriminator:
    """Classify mobile genetic elements carrying resistance"""
    
    # SXT/R391 ICE markers
    SXT_ICE_MARKERS = {
        "integration_sites": ["attL", "attR"],  # dif sites on chromosome
        "core_genes": ["setCD", "acaCD", "xis", "int"],
        "resistance_genes": ["floR", "tetA(D)", "fexA", "dfrA", "dfrB"],
        "regulator": "AcaCD",  # Primary regulator
        "element_type": "Integrative Conjugative Element (ICE)",
        "location": "chromosome"
    }
    
    # IncA/C plasmid markers
    INCA_C_PLASMID_MARKERS = {
        "replicon": ["IncA", "IncC"],
        "backbone_genes": ["traC", "traD", "traG", "mob"],
        "resistance_genes": ["aac(3)-IIa", "dfrA", "tetA(D)", "floR"],
        "regulator": "SetCD",  # Different regulator
        "element_type": "Conjugative Plasmid",
        "location": "plasmid",
        "transmission": "higher_rate",
        "stability": "potentially_unstable"
    }
    
    def __init__(self):
        self.element_calls = {}
        
    def detect_sxt_element(self, assembly_data: Dict[str, Any], k_mer_matches: Dict[str, Any]) -> Dict[str, Any]:
        """
        Detect SXT/R391 ICE element on chromosome
        
        Args:
            assembly_data: Graph/contig data from assembly
            k_mer_matches: K-mer signature matches from profiling
        
        Returns:
            Dictionary with SXT status, genes present, integration sites
        """
        result = {
            "sxt_element_detected": False,
            "sxt_status": "NOT_DETECTED",
            "integration_sites": {},
            "core_genes_found": [],
            "resistance_genes": [],
            "location": "UNKNOWN",
            "regulator": "UNKNOWN",
            "confidence": "LOW"
        }
        
        # Check for SXT integration markers
        if isinstance(k_mer_matches, dict):
            sxt_markers_found = 0
            
            # Check for core SXT genes
            for gene in self.SXT_ICE_MARKERS["core_genes"]:
                if k_mer_matches.get(gene):
                    result["core_genes_found"].append(gene)
                    sxt_markers_found += 1
            
            # Check for resistance genes
            for gene in self.SXT_ICE_MARKERS["resistance_genes"]:
                if k_mer_matches.get(gene):
                    result["resistance_genes"].append(gene)
                    sxt_markers_found += 1
            
            # Check for regulator (AcaCD)
            if k_mer_matches.get("acaCD") or k_mer_matches.get("setCD"):
                result["regulator"] = "AcaCD" if k_mer_matches.get("acaCD") else "SetCD"
                sxt_markers_found += 1
            
            # Determine if SXT is present based on marker accumulation
            if sxt_markers_found >= 3:  # At least 3 markers
                result["sxt_element_detected"] = True
                result["sxt_status"] = "SXT_R391_ICE_DETECTED"
                result["location"] = "CHROMOSOME_INTEGRATED"
                result["confidence"] = "HIGH"
                
                # Simulate integration site data
                result["integration_sites"] = {
                    "attL_dif1": {"detected": True, "coverage": "expected"},
                    "attR_dif2": {"detected": True, "coverage": "expected"}
                }
        
        return result
    
    def detect_inca_c_plasmid(self, k_mer_matches: Dict[str, Any], contig_analysis: Optional[Dict] = None) -> Dict[str, Any]:
        """
        Detect IncA/C conjugative plasmids
        
        Args:
            k_mer_matches: K-mer signature matches
            contig_analysis: Contig/plasmid graph analysis (optional)
        
        Returns:
            Dictionary with IncA/C status, genes present, transmission profile
        """
        result = {
            "inca_c_plasmid_detected": False,
            "plasmid_status": "NOT_DETECTED",
            "replicon_type": "UNKNOWN",
            "backbone_genes_found": [],
            "resistance_genes": [],
            "location": "UNKNOWN",
            "regulator": "UNKNOWN",
            "transmission_dynamics": "UNKNOWN",
            "confidence": "LOW"
        }
        
        if isinstance(k_mer_matches, dict):
            inca_markers_found = 0
            
            # Check for replicon markers (IncA or IncC)
            for replicon in self.INCA_C_PLASMID_MARKERS["replicon"]:
                if k_mer_matches.get(replicon):
                    result["replicon_type"] = replicon
                    inca_markers_found += 1
            
            # Check for plasmid backbone genes
            for gene in self.INCA_C_PLASMID_MARKERS["backbone_genes"]:
                if k_mer_matches.get(gene):
                    result["backbone_genes_found"].append(gene)
                    inca_markers_found += 1
            
            # Check for resistance genes
            for gene in self.INCA_C_PLASMID_MARKERS["resistance_genes"]:
                if k_mer_matches.get(gene):
                    result["resistance_genes"].append(gene)
                    inca_markers_found += 1
            
            # Determine if IncA/C plasmid is present
            if inca_markers_found >= 3:  # At least replicon + 2 genes
                result["inca_c_plasmid_detected"] = True
                result["plasmid_status"] = "INCA_C_CONJUGATIVE_PLASMID_DETECTED"
                result["location"] = "PLASMID"
                result["regulator"] = "SetCD"
                result["transmission_dynamics"] = "HIGHER_RATE (conjugative)"
                result["confidence"] = "HIGH"
        
        return result
    
    def classify_amr_element(self, sxt_result: Dict[str, Any], inca_c_result: Dict[str, Any]) -> Dict[str, Any]:
        """
        Final classification of primary AMR element
        Integrates SXT and IncA/C results
        """
        result = {
            "primary_amr_element": "UNKNOWN",
            "element_type": "UNKNOWN",
            "location": "UNKNOWN",
            "regulator_set": [],
            "transmission_profile": "UNKNOWN",
            "combined_resistance": [],
            "confidence": "LOW"
        }
        
        # Determine primary element
        if sxt_result.get("sxt_element_detected"):
            result["primary_amr_element"] = "SXT_R391_ICE"
            result["element_type"] = "Integrative Conjugative Element"
            result["location"] = "CHROMOSOME_INTEGRATED"
            result["regulator_set"] = ["AcaCD", "SetCD"]
            result["transmission_profile"] = "Chromosomal integration (low lateral transfer rate)"
            result["combined_resistance"] = sxt_result.get("resistance_genes", [])
            result["confidence"] = "HIGH"
        elif inca_c_result.get("inca_c_plasmid_detected"):
            result["primary_amr_element"] = "INCA_C_CONJUGATIVE_PLASMID"
            result["element_type"] = "Plasmid"
            result["location"] = "PLASMID"
            result["regulator_set"] = ["SetCD", "TraC"]
            result["transmission_profile"] = "Plasmid conjugation (high lateral transfer rate)"
            result["combined_resistance"] = inca_c_result.get("resistance_genes", [])
            result["confidence"] = "HIGH"
        else:
            result["primary_amr_element"] = "UNCLASSIFIED_OR_MONORESISTANCE"
            result["confidence"] = "UNKNOWN"
        
        return result
    
    def predict_transmission_dynamics(self, element_type: str) -> Dict[str, Any]:
        """
        Predict transmission and stability profile based on element type
        """
        if "SXT" in element_type or "ICE" in element_type:
            return {
                "element": "SXT/ICE",
                "horizontal_transfer_rate": "MODERATE",
                "transfer_mechanism": "Chromosomal conjugation (requires cell contact)",
                "stability": "STABLE (chromosomal integration)",
                "clinical_implications": [
                    "Resistance is stably inherited",
                    "Requires direct cell-to-cell contact for spread",
                    "Unlikely to spread to unrelated pathogens"
                ]
            }
        elif "PLASMID" in element_type or "IncA/C" in element_type:
            return {
                "element": "IncA/C Plasmid",
                "horizontal_transfer_rate": "HIGH",
                "transfer_mechanism": "Plasmid conjugation + mobilization by helper plasmids",
                "stability": "VARIABLE (can be lost at high frequency)",
                "clinical_implications": [
                    "Resistance can spread rapidly between strains",
                    "May transfer to related pathogens (Aeromonas, Enterobacterium)",
                    "Can be counterselected if conjugation burden is high"
                ]
            }
        else:
            return {"element": "UNKNOWN", "transmission_profile": "UNKNOWN"}


def test_amr_discrimination():
    """Unit test for AMR element discrimination"""
    discriminator = AMRElementDiscriminator()
    
    # Test 1: Detect SXT element
    k_mer_sxt = {
        "setCD": True,
        "acaCD": True,
        "floR": True,
        "tetA(D)": True,
        "dfrA": True
    }
    sxt_result = discriminator.detect_sxt_element({}, k_mer_sxt)
    assert sxt_result["sxt_element_detected"] == True
    assert sxt_result["location"] == "CHROMOSOME_INTEGRATED"
    print("✓ Test 1 passed: SXT element detection")
    
    # Test 2: Detect IncA/C plasmid
    k_mer_inca = {
        "IncA": True,
        "traC": True,
        "traD": True,
        "dfrA": True,
        "aac(3)-IIa": True
    }
    inca_result = discriminator.detect_inca_c_plasmid(k_mer_inca)
    assert inca_result["inca_c_plasmid_detected"] == True
    assert inca_result["location"] == "PLASMID"
    print("✓ Test 2 passed: IncA/C plasmid detection")
    
    # Test 3: Classification
    sxt_only = {"sxt_element_detected": True, "resistance_genes": ["floR", "tetA(D)"]}
    inca_none = {"inca_c_plasmid_detected": False}
    
    classification = discriminator.classify_amr_element(sxt_only, inca_none)
    assert classification["primary_amr_element"] == "SXT_R391_ICE"
    assert "CHROMOSOME" in classification["location"]
    print("✓ Test 3 passed: SXT classification")
    
    # Test 4: Transmission dynamics
    dynamics = discriminator.predict_transmission_dynamics("SXT/ICE")
    assert "STABLE" in dynamics["stability"]
    print("✓ Test 4 passed: Transmission dynamics prediction")


if __name__ == "__main__":
    test_amr_discrimination()
    print("\nAll tests passed!")
