#!/usr/bin/env python3
"""
Dual-Marker Serogroup Detection (O1 vs O139 vs Phage Scar)

Logic:
- O1 El Tor: rfb present, wbf absent
- O139 Bengal: rfb absent, wbf present (22kb deletion + 35kb insertion)
- Phage Scar: rfb absent, wbf absent (predation damage, not O139)
"""

import argparse
import subprocess
import json
import sys
from pathlib import Path


def get_region_coverage(bam_path: str, region: str) -> float:
    """
    Calculate average coverage for a genomic region.
    
    Args:
        bam_path: Path to sorted BAM file
        region: Region string (e.g., "CP003069.1:100-200")
    
    Returns:
        Average depth across region
    """
    cmd = ["samtools", "depth", "-a", "-r", region, bam_path]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        depths = []
        
        for line in result.stdout.splitlines():
            parts = line.split('\t')
            if len(parts) >= 3:
                depths.append(int(parts[2]))
        
        if not depths:
            return 0.0
        
        return sum(depths) / len(depths)
    
    except subprocess.CalledProcessError as e:
        print(f"Error calculating coverage: {e}", file=sys.stderr)
        return 0.0


def check_bridging_reads(bam_path: str, region_start: str, region_end: str) -> bool:
    """
    Check if there are reads that bridge a gap (span both flanks).
    This distinguishes real deletions from coverage dropouts.
    
    Args:
        bam_path: Path to sorted BAM file
        region_start: Upstream flank region
        region_end: Downstream flank region
    
    Returns:
        True if bridging reads exist
    """
    # Extract read names from both flanks
    cmd_start = f"samtools view {bam_path} {region_start} | cut -f1 | sort -u"
    cmd_end = f"samtools view {bam_path} {region_end} | cut -f1 | sort -u"
    
    try:
        reads_start = set(subprocess.check_output(cmd_start, shell=True, text=True).splitlines())
        reads_end = set(subprocess.check_output(cmd_end, shell=True, text=True).splitlines())
        
        # Bridging reads appear in both sets
        bridging = reads_start & reads_end
        return len(bridging) > 0
    
    except subprocess.CalledProcessError:
        return False


def call_wbeT_variants(bam_path: str, reference: str, region: str) -> list:
    """
    Call variants in wbeT using bcftools to detect Inaba switch (Stop Codons/SNPs).
    """
    try:
        # mpileup -> call -> query
        # We look for ANY non-reference homozygous or heterozygous call in this gene
        cmd = (
            f"bcftools mpileup -Ou -f {reference} -r {region} {bam_path} | "
            f"bcftools call -mv -Ou | "
            "bcftools query -f '%POS:%REF>%ALT\\n'"
        )
        output = subprocess.check_output(cmd, shell=True, text=True).strip()
        if not output:
            return []
        
        return output.splitlines()
    except subprocess.CalledProcessError:
        return []


def detect_serogroup(bam_path: str, reference: str) -> dict:
    """
    Detect serogroup using dual-marker logic.
    
    Regions (Haiti 2010EL-1786 coordinates):
    - rfb cluster: CP003069.1:370000-390000 (O1 antigen)
    - wbf cluster: Not in O1 reference (must check against O139 ref)
    - Flanking genes: gmhD, rjg (for deletion confirmation)
    """
    # Define marker regions based on validated coordinates (Haiti 2010EL-1786)
    # The 'Inaba Trap': Inaba has a mutated/pseudogenized wbeT but IS O1.
    # We must separate Serogroup (O1) from Serotype (Ogawa/Inaba).
    
    # 1. Serogroup Anchors (Flanking Genes) - Must be present for any O1
    gmhD_region = "CP003069.1:369000-370000"   # Upstream anchor
    rjg_region = "CP003069.1:390000-391000"    # Downstream anchor
    
    # 2. Serotype Marker (wbeT/rfbT) - The differentiator
    # Inaba often has mutations here, but coverage should still exist unless deleted.
    # If absent/low coverage but flanks exist -> Likely Inaba or partial deletion.
    wbeT_region = "CP003069.1:2678186-2678980" # Using rfbH/wbeT candidate region
    
    # O139 check remains similar (bridging reads or wbf)
    rfb_region = "CP003069.1:370000-390000" # Full O1 cluster
    
    # Calculate coverage
    gmhD_cov = get_region_coverage(bam_path, gmhD_region)
    rjg_cov = get_region_coverage(bam_path, rjg_region)
    wbeT_cov = get_region_coverage(bam_path, wbeT_region)
    rfb_cluster_cov = get_region_coverage(bam_path, rfb_region)
    
    # Check for bridging reads (confirms deletion for O139)
    # Bridging reads span from gmhD to rjg, skipping the O1 cluster
    has_bridging = check_bridging_reads(bam_path, gmhD_region, rjg_region)
    
    # Logic Decision Tree
    serogroup = "UNKNOWN"
    serotype = "UNKNOWN"
    confidence = "LOW"
    reason = []
    
    # Step 1: Serogroup Determination via Flanks
    is_o1_backbone = (gmhD_cov > 5) and (rjg_cov > 5)
    
    if is_o1_backbone:
        if rfb_cluster_cov > 5:
            serogroup = "O1"
            confidence = "HIGH"
            reason.append(f"O1 Backbone confirmed (gmhD={gmhD_cov:.1f}x, rjg={rjg_cov:.1f}x).")
            reason.append(f"O1 Antigen cluster present ({rfb_cluster_cov:.1f}x).")
            
            # Step 2: Serotype Determination via wbeT
            # Ogawa = Functional wbeT. Inaba = Mutated/Non-functional wbeT.
            # In sequencing, 'inaba' often looks like 'variant wbeT' or 'low mapping' if divergent.
            # If wbeT coverage is good, it's likely Ogawa (or Inaba with SNP). 
            # If wbeT is missing but rest of cluster is there -> Strong Inaba candidate (or rare Hikojima).
            
            if wbeT_cov > 5:
                # Genomic wbeT is present. Check for SNPs (Premature Stop / Inaba Switch).
                # Specifically checking for mutations in wbeT region.
                variants = call_wbeT_variants(bam_path, reference, wbeT_region)
                if variants:
                    serotype = "Inaba (Confirmed SNP)"
                    reason.append(f"wbeT gene present ({wbeT_cov:.1f}x) but contains mutations: {', '.join(variants)}.")
                else:
                    serotype = "Ogawa (Wild Type)"
                    reason.append(f"wbeT gene present ({wbeT_cov:.1f}x) and Wild Type (Ogawa).")
            else:
                # wbeT maps poorly or is deleted, but cluster is there.
                serotype = "Inaba (Likely)"
                reason.append(f"wbeT gene absent/divergent ({wbeT_cov:.1f}x) despite O1 cluster presence.")
                
        elif has_bridging:
            # Flanks present, cluster missing, bridging reads found -> O139 (or O139-like replacement)
            serogroup = "O139"
            serotype = "Bengal"
            confidence = "HIGH"
            reason.append("O1 Backbone present but Antigen cluster replaced (Bridging Reads found).")
            
        else:
            # Flanks present, cluster missing, no bridging.
            serogroup = "NON_O1_NON_O139" # Possible Phage Scar or rough strain
            serotype = "Rough/Phage-Damaged"
            confidence = "MEDIUM"
            reason.append(f"O1 Backbone present but Antigen cluster missing ({rfb_cluster_cov:.1f}x). No bridging reads.")
            
    else:
        # No O1 backbone
        serogroup = "NON_O1"
        confidence = "HIGH"
        reason.append(f"O1 Backbone absent (gmhD={gmhD_cov:.1f}x, rjg={rjg_cov:.1f}x).")
    
    return {
        "serogroup": serogroup,
        "serotype": serotype,
        "confidence": confidence,
        "reason": " ".join(reason),
        "metrics": {
            "gmhD_coverage": round(gmhD_cov, 2),
            "rjg_coverage": round(rjg_cov, 2),
            "rfb_cluster_coverage": round(rfb_cluster_cov, 2),
            "wbeT_coverage": round(wbeT_cov, 2),
            "has_bridging_reads": has_bridging
        }
    }


def main():
    parser = argparse.ArgumentParser(description="Detect Vibrio cholerae serogroup")
    parser.add_argument("--bam", required=True, help="Sorted BAM file")
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument("--output", required=True, help="Output JSON report")
    
    args = parser.parse_args()
    
    print(f"🧬 Detecting serogroup from {Path(args.bam).name}...")
    
    result = detect_serogroup(args.bam, args.reference)
    
    # Save report
    with open(args.output, 'w') as f:
        json.dump(result, f, indent=2)
    
    print(f"   ✅ Serogroup: {result['serogroup']} ({result['confidence']} confidence)")
    print(f"   📊 {result['reason']}")
    print(f"   📁 Report saved: {args.output}")


if __name__ == "__main__":
    main()
