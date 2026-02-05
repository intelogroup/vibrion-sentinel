#!/usr/bin/env python3
"""
Pipeline Mode Selector: Field Rapid vs Laboratory Full

Optimizes consensus generation for resource-constrained field deployment
vs comprehensive laboratory analysis.
"""

from typing import Dict


def configure_pipeline_depth(mode: str) -> Dict:
    """
    Configure pipeline parameters based on deployment mode.
    
    Args:
        mode: "FIELD_RAPID" or "LABORATORY_FULL"
    
    Returns:
        Configuration dictionary
    """
    if mode == "FIELD_RAPID":
        return {
            "polish_rounds": 1,
            "local_assembly": False,
            "sxt_assembly": False,  # Skip SPAdes in field mode
            "maf_threshold": None,  # Skip minority variant tracking
            "min_depth": 3,  # Lower threshold for masking
            "checksum_validation": False,
            "estimated_runtime_min": 15,
            "impact": "Generate Rapid Alert only",
            "description": "Optimized for field deployment with limited compute"
        }
    
    elif mode == "LABORATORY_FULL":
        return {
            "polish_rounds": 2,
            "local_assembly": True,  # SPAdes for SXT region
            "sxt_assembly": True,  # Full de novo SXT assembly
            "maf_threshold": 0.20,  # Track minority alleles ≥20%
            "min_depth": 5,  # Strict masking threshold
            "checksum_validation": True,  # Validate housekeeping genes
            "estimated_runtime_min": 90,
            "impact": "Full Forensic Publication-Ready Genome",
            "description": "Comprehensive analysis with all validation layers"
        }
    
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'FIELD_RAPID' or 'LABORATORY_FULL'")


def get_mode_from_config(config_path: str = None) -> str:
    """
    Read mode from config file or environment.
    Defaults to LABORATORY_FULL if not specified.
    """
    # TODO: Read from config.yaml or environment variable
    # For now, default to full mode
    return "LABORATORY_FULL"


if __name__ == "__main__":
    import json
    import sys
    
    mode = sys.argv[1] if len(sys.argv) > 1 else "LABORATORY_FULL"
    
    config = configure_pipeline_depth(mode)
    print(json.dumps(config, indent=2))
