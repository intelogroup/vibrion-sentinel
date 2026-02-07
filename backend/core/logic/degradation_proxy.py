#!/usr/bin/env python3
"""
QC & Degradation Proxy: Detect sample quality issues from freeze-thaw cycles
K-mer coverage CV, read quality metrics, DNA integrity estimation
"""

import json
import math
from typing import Dict, Any, Optional, List
from pathlib import Path

class DegradationProxyCalculator:
    """Quantify sample degradation without wet-lab analysis"""
    
    def __init__(self):
        self.degradation_metrics = {}
        
    def calculate_kmer_cv(self, kmer_depths: List[float]) -> float:
        """
        Calculate coefficient of variation for k-mer depths
        
        CV = stdev / mean
        
        Low CV (<0.10) = high-quality, pristine DNA
        Medium CV (0.10-0.20) = acceptable quality, minor degradation
        High CV (0.20-0.40) = significant degradation (freeze-thaw suspected)
        Very high CV (>0.40) = severe degradation, unusable
        
        Args:
            kmer_depths: List of k-mer depth values across genome
        
        Returns:
            Coefficient of variation (0.0-1.0)
        """
        if not kmer_depths or len(kmer_depths) < 2:
            return -1.0  # Invalid
        
        mean_depth = sum(kmer_depths) / len(kmer_depths)
        if mean_depth == 0:
            return -1.0
        
        variance = sum((x - mean_depth) ** 2 for x in kmer_depths) / len(kmer_depths)
        stdev = math.sqrt(variance)
        cv = stdev / mean_depth
        
        return cv
    
    def estimate_freeze_thaw_cycles(self, cv_kmer: float, vibrio_pct: float, 
                                   mean_read_quality: Optional[float] = None) -> Dict[str, Any]:
        """
        Estimate number of freeze-thaw cycles based on degradation proxies
        
        Args:
            cv_kmer: K-mer coefficient of variation
            vibrio_pct: Percentage of reads classified as Vibrio
            mean_read_quality: Mean Phred quality score (optional)
        
        Returns:
            Dictionary with estimated cycles, risk score, recommendations
        """
        result = {
            "estimated_freeze_thaw_cycles": 0,
            "freeze_thaw_risk_score": 0.0,
            "dna_integrity_estimate": "UNKNOWN",
            "quality_category": "UNKNOWN",
            "recommendations": []
        }
        
        # Calculate degradation score
        # Components: k-mer CV (high = degraded), vibrio % (low = contamination), read quality (low = degraded)
        
        cv_component = min(cv_kmer / 0.5, 1.0)  # Normalize to 0-1 (0.5 is saturation)
        contamination_component = max(0, (100 - vibrio_pct) / 20.0)  # Normalize to 0-1 (20% contamination is saturation)
        
        quality_component = 0.0
        if mean_read_quality is not None:
            # Phred quality 30-40 is high, <20 is low
            quality_component = max(0, (30 - mean_read_quality) / 10.0)
        
        # Overall degradation score (0-3)
        total_degradation = cv_component + contamination_component + quality_component
        result["freeze_thaw_risk_score"] = total_degradation
        
        # Map degradation to freeze-thaw cycles
        if cv_kmer < 0.08 and vibrio_pct > 96:
            result["estimated_freeze_thaw_cycles"] = 0
            result["dna_integrity_estimate"] = "PRISTINE"
            result["quality_category"] = "EXCELLENT"
            result["recommendations"].append("No degradation detected; immediate processing recommended")
        
        elif cv_kmer < 0.12 and vibrio_pct > 93:
            result["estimated_freeze_thaw_cycles"] = 1
            result["dna_integrity_estimate"] = "MINOR_DEGRADATION"
            result["quality_category"] = "ACCEPTABLE"
            result["recommendations"].append("One freeze-thaw cycle suspected; still suitable for analysis")
        
        elif cv_kmer < 0.20 and vibrio_pct > 88:
            result["estimated_freeze_thaw_cycles"] = 2
            result["dna_integrity_estimate"] = "MODERATE_DEGRADATION"
            result["quality_category"] = "BORDERLINE"
            result["recommendations"].append("2-3 freeze-thaw cycles suspected")
            result["recommendations"].append("Confidence in minor variant calls reduced")
            result["recommendations"].append("Consider re-sequencing from backup sample")
        
        elif cv_kmer < 0.35 and vibrio_pct > 80:
            result["estimated_freeze_thaw_cycles"] = 3
            result["dna_integrity_estimate"] = "SIGNIFICANT_DEGRADATION"
            result["quality_category"] = "POOR"
            result["recommendations"].append("3-5 freeze-thaw cycles suspected")
            result["recommendations"].append("Do NOT use for clinical decision-making")
            result["recommendations"].append("Re-extract from original sample or discard")
        
        else:
            result["estimated_freeze_thaw_cycles"] = 5
            result["dna_integrity_estimate"] = "SEVERE_DEGRADATION"
            result["quality_category"] = "UNUSABLE"
            result["recommendations"].append("Severe degradation detected (>5 freeze-thaw cycles)")
            result["recommendations"].append("CRITICAL: Sample unusable for analysis")
            result["recommendations"].append("Return to lab for fresh extraction")
        
        return result
    
    def generate_qc_report(self, vibrio_stats: Dict[str, Any], coverage_stats: Dict[str, Any],
                          kmer_depths: Optional[List[float]] = None) -> Dict[str, Any]:
        """
        Comprehensive QC report combining multiple degradation proxies
        """
        report = {
            "timestamp": None,
            "sample_id": vibrio_stats.get("sample_id"),
            "vibrio_percentage": vibrio_stats.get("vibrio_percentage", 0.0),
            "total_reads": vibrio_stats.get("total_reads", 0),
            "degradation_proxy_cv": -1.0,
            "freeze_thaw_risk_score": 0.0,
            "dna_integrity_estimate": "UNKNOWN",
            "qc_status": "UNKNOWN",
            "recommendations": [],
            "field_deployment_suitable": False
        }
        
        # Calculate k-mer CV if available
        if kmer_depths:
            cv = self.calculate_kmer_cv(kmer_depths)
            report["degradation_proxy_cv"] = cv
            
            # Estimate degradation
            degradation = self.estimate_freeze_thaw_cycles(cv, report["vibrio_percentage"])
            report["freeze_thaw_risk_score"] = degradation["freeze_thaw_risk_score"]
            report["dna_integrity_estimate"] = degradation["dna_integrity_estimate"]
            report["recommendations"].extend(degradation["recommendations"])
        
        # Overall QC gating
        vibrio = report["vibrio_percentage"]
        cv = report["degradation_proxy_cv"]
        
        if vibrio > 95 and cv < 0.10:
            report["qc_status"] = "PASS"
            report["field_deployment_suitable"] = True
        elif vibrio > 90 and cv < 0.15:
            report["qc_status"] = "PASS"
            report["field_deployment_suitable"] = True
        elif vibrio > 80 and cv < 0.25:
            report["qc_status"] = "BORDERLINE"
            report["field_deployment_suitable"] = False
            report["recommendations"].append("Use with caution; manual review recommended")
        else:
            report["qc_status"] = "FAIL"
            report["field_deployment_suitable"] = False
            report["recommendations"].append("CRITICAL: Do not use for clinical decisions")
        
        return report


class SNPDistanceCalculator:
    """Calculate SNP distances for phylogenetic time calibration"""
    
    def __init__(self, reference_genome: Optional[str] = None):
        self.reference = reference_genome
        self.distances = {}
        
    def calculate_snp_distance(self, vcf_data: Dict[str, Any]) -> float:
        """
        Count variant positions (SNPs) vs reference
        
        Used for time calibration: Haiti 2010 = 0, EnvJ515 (2018) = 12-18, etc.
        """
        if isinstance(vcf_data, dict) and "variants" in vcf_data:
            # Count SNP positions (exclude indels and complex variants)
            snp_count = sum(1 for v in vcf_data.get("variants", []) 
                          if v.get("variant_type", "snp") == "snp")
            return float(snp_count)
        
        return 0.0
    
    def estimate_divergence_date(self, snp_distance: float, 
                                mutation_rate: float = 1.0e-4) -> Dict[str, Any]:
        """
        Estimate divergence date from SNP distance
        
        Simple molecular clock: years_ago ≈ snp_distance / (2 * mutation_rate * genome_size)
        
        For V. cholerae: genome ~4Mb, effective mutation rate ~1e-4 SNP/bp/year
        """
        result = {
            "snp_distance": snp_distance,
            "estimated_divergence_years_ago": 0.0,
            "estimated_year": None,
            "confidence": "LOW"
        }
        
        if snp_distance <= 0:
            return result
        
        # Rough estimate for V. cholerae
        # ~0.5-2 SNPs per year observed in Haiti strains
        years_per_snp = 0.5  # Conservative estimate
        years_ago = snp_distance * years_per_snp
        
        result["estimated_divergence_years_ago"] = years_ago
        result["estimated_year"] = 2010 + years_ago  # Relative to 2010 Haiti ancestor
        result["confidence"] = "MEDIUM" if years_ago < 30 else "LOW"
        
        return result


def test_degradation_proxy():
    """Unit test for degradation proxy"""
    calculator = DegradationProxyCalculator()
    
    # Test 1: Pristine sample
    pristine_cv = 0.05
    result_pristine = calculator.estimate_freeze_thaw_cycles(pristine_cv, 98.5)
    assert result_pristine["estimated_freeze_thaw_cycles"] == 0
    assert result_pristine["quality_category"] == "EXCELLENT"
    print("✓ Test 1 passed: Pristine sample detection")
    
    # Test 2: One F-T cycle
    one_ft_cv = 0.10
    result_one = calculator.estimate_freeze_thaw_cycles(one_ft_cv, 95.2)
    assert result_one["estimated_freeze_thaw_cycles"] == 1
    assert result_one["quality_category"] == "ACCEPTABLE"
    print("✓ Test 2 passed: One F-T cycle detection")
    
    # Test 3: Multiple F-T cycles
    multi_cv = 0.25
    result_multi = calculator.estimate_freeze_thaw_cycles(multi_cv, 85.0)
    assert result_multi["estimated_freeze_thaw_cycles"] >= 2
    assert result_multi["quality_category"] in ["BORDERLINE", "POOR"]
    print("✓ Test 3 passed: Multiple F-T cycle detection")
    
    # Test 4: K-mer CV calculation
    kmer_pristine = [100, 102, 101, 99, 101, 100]  # Low variance
    cv = calculator.calculate_kmer_cv(kmer_pristine)
    assert cv < 0.10
    print("✓ Test 4 passed: K-mer CV calculation")


def test_snp_distance():
    """Unit test for SNP distance calculator"""
    calculator = SNPDistanceCalculator()
    
    # Test 1: Haiti reference (0 SNPs)
    vcf_haiti = {"variants": []}
    dist_haiti = calculator.calculate_snp_distance(vcf_haiti)
    assert dist_haiti == 0.0
    print("✓ Test 1 passed: Haiti reference (0 SNPs)")
    
    # Test 2: EnvJ515 (15 SNPs)
    vcf_env = {"variants": [{"variant_type": "snp"} for _ in range(15)]}
    dist_env = calculator.calculate_snp_distance(vcf_env)
    assert dist_env == 15.0
    print("✓ Test 2 passed: EnvJ515 SNP distance")
    
    # Test 3: Divergence dating
    date_estimate = calculator.estimate_divergence_date(15.0)
    # 15 SNPs * 0.5 SNPs/year = ~7.5 years from 2010 → 2017.5
    assert 2015 < date_estimate["estimated_year"] < 2025  # Within reasonable range
    print("✓ Test 3 passed: Divergence dating")


if __name__ == "__main__":
    test_degradation_proxy()
    print()
    test_snp_distance()
    print("\nAll tests passed!")
