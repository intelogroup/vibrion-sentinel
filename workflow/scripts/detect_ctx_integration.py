#!/usr/bin/env python3
"""
CTXφ Dual-Site Integration Detector
Detects cholera toxin prophage integration at Chr1 dif1 and Chr2 dif2 sites
"""

import argparse
import json
import subprocess
import sys
import os
import tempfile
from pathlib import Path
from Bio import SeqIO
import pysam

# CTX integration site coordinates (2010EL-1786 reference)
DEFAULT_DIF_SITES = {
    "dif1_chr1": {
        "chrom": "CP003069.1",
        "start": 1041000,
        "end": 1047000,
        "description": "Chromosome 1 CTXphi island (canonical Haiti location)",
        "query": "ctxB"
    },
    "dif2_chr2": {
        "chrom": "CP003070.1",
        "start": 880000,
        "end": 885000,
        "description": "Chromosome 2 pspE region (atypical integration site)",
        "query": "pspE" # use pspE as proxy for Chr2 site
    }
}

def discover_region_via_blast(query_name, target_fasta):
    """
    Search for the locus using BLAST to find coordinates in the current reference.
    """
    ref_loci_path = "data/references/reference_loci.fasta"
    if not os.path.exists(ref_loci_path):
        return None

    # 1. Extract the query sequence from reference_loci.fasta
    query_seq = None
    for record in SeqIO.parse(ref_loci_path, "fasta"):
        if record.id == query_name:
            query_seq = str(record.seq)
            break
    
    if not query_seq:
        return None

    # 2. BLAST against target
    with tempfile.NamedTemporaryFile(suffix=".fasta", mode="w") as tmp_query:
        tmp_query.write(f">{query_name}\n{query_seq}\n")
        tmp_query.flush()

        cmd = [
            "blastn",
            "-query", tmp_query.name,
            "-subject", str(target_fasta),
            "-outfmt", "6 sseqid sstart send pident length",
            "-perc_identity", "70",
            "-word_size", "7",
            "-max_target_seqs", "1"
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            hits = result.stdout.strip().split("\n")
            if not hits or not hits[0]:
                return None
            
            # 3. Extract the hit coordinates
            parts = hits[0].split("\t")
            sseqid, sstart, send = parts[0], int(parts[1]), int(parts[2])
            
            # Ensure sstart < send
            start, end = min(sstart, send), max(sstart, send)
            # Expand a bit to cover the integration site
            return {"chrom": sseqid, "start": max(0, start - 2000), "end": end + 2000}
        except Exception:
            return None

def calculate_depth(bam_file: Path, chrom: str, start: int, end: int) -> dict:
    """
    Calculate mean depth and detect integration spike.
    """
    # Use samtools depth for the region
    cmd = [
        "samtools", "depth",
        "-r", f"{chrom}:{start}-{end}",
        str(bam_file)
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        depths = []
        
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                depths.append(int(parts[2]))
        
        if not depths:
            return {"mean_depth": 0, "max_depth": 0, "spike_detected": False}
        
        mean_depth = sum(depths) / len(depths)
        max_depth = max(depths)
        
        # Calculate flanking region depth for comparison
        flank_start = max(0, start - 5000)
        flank_end = max(0, start - 1000)
        flank_cmd = [
            "samtools", "depth",
            "-r", f"{chrom}:{flank_start}-{flank_end}",
            str(bam_file)
        ]
        flank_result = subprocess.run(flank_cmd, capture_output=True, text=True, check=True)
        flank_depths = [int(line.split('\t')[2]) for line in flank_result.stdout.strip().split('\n') if line]
        flank_mean = sum(flank_depths) / len(flank_depths) if flank_depths else mean_depth
        
        # Spike detected if region depth > 3x flanking
        spike_detected = mean_depth > (flank_mean * 3) if flank_mean > 0 else False
        
        return {
            "mean_depth": round(mean_depth, 2),
            "max_depth": max_depth,
            "spike_detected": spike_detected,
            "flanking_depth": round(flank_mean, 2)
        }
    
    except subprocess.CalledProcessError:
        print(f"Warning: Could not calculate depth for {chrom}:{start}-{end}", file=sys.stderr)
        return {"mean_depth": 0, "max_depth": 0, "spike_detected": False}

def count_bridging_reads(bam_file: Path, chrom: str, start: int, end: int) -> int:
    """
    Count reads that span the CTX junction (bridging reads).
    """
    cmd = [
        "samtools", "view",
        "-c",  # count only
        "-F", "4",  # exclude unmapped
        str(bam_file),
        f"{chrom}:{start}-{end}"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return int(result.stdout.strip())
    except (subprocess.CalledProcessError, ValueError):
        return 0

def detect_ctx_integration(bam_file: Path, dif_sites, min_depth: int = 5) -> dict:
    """
    Detect CTX prophage integration at both dif sites.
    """
    results = {}
    total_copies = 0
    
    for site_name, coords in dif_sites.items():
        if coords is None:
            results[site_name] = {"detected": False, "mean_depth": 0, "copy_estimate": 0, "bridging_reads": 0, "spike_detected": False}
            continue

        print(f"   Checking {site_name}...")
        
        depth_data = calculate_depth(
            bam_file,
            coords["chrom"],
            coords["start"],
            coords["end"]
        )
        
        bridging_reads = count_bridging_reads(
            bam_file,
            coords["chrom"],
            coords["start"],
            coords["end"]
        )
        
        # Integration detected if:
        # 1. Mean depth >= min_depth
        # 2. Depth spike detected OR bridging reads present
        detected = (
            depth_data["mean_depth"] >= min_depth and
            (depth_data["spike_detected"] or bridging_reads > 10)
        )
        
        # Estimate copy number based on depth ratio
        copy_estimate = 0
        if detected:
            copy_estimate = 1
            if depth_data["mean_depth"] > depth_data.get("flanking_depth", 50) * 5:
                copy_estimate = 2  # Possible tandem copies
        
        results[site_name] = {
            "detected": detected,
            "mean_depth": depth_data["mean_depth"],
            "copy_estimate": copy_estimate,
            "bridging_reads": bridging_reads,
            "spike_detected": depth_data["spike_detected"],
            "region": f"{coords['chrom']}:{coords['start']}-{coords['end']}"
        }
        
        total_copies += copy_estimate
    
    # Determine overall CTX status
    if total_copies == 0:
        ctx_status = "NOT_DETECTED"
        warning = None
    elif results["dif1_chr1"]["detected"] and results["dif2_chr2"]["detected"]:
        ctx_status = "MULTIPLE_SITES_DETECTED"
        warning = "⚠️ CTX prophage detected at BOTH Chr1 and Chr2 dif sites. Atypical El Tor strain with potential for higher toxin production."
    elif results["dif1_chr1"]["detected"]:
        ctx_status = "INTEGRATED_CHR1"
        warning = None
    elif results["dif2_chr2"]["detected"]:
        ctx_status = "INTEGRATED_CHR2"
        warning = "⚠️ CTX prophage at Chr2 only. Unusual integration pattern."
    else:
        ctx_status = "AMBIGUOUS"
        warning = "Low coverage or ambiguous integration signal"
    
    return {
        "ctx_status": ctx_status,
        "integration_sites": results,
        "total_copy_number": total_copies,
        "warning": warning
    }

def main():
    parser = argparse.ArgumentParser(
        description="Detect CTX prophage integration at Chr1/Chr2 dif sites"
    )
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument("--output", required=True, help="Output JSON report")
    parser.add_argument("--sample", help="Sample ID")
    parser.add_argument("--min-depth", type=int, default=5,
                       help="Minimum depth for integration call (default: 5)")
    
    args = parser.parse_args()
    
    print("🔬 CTXφ Dual-Site Integration Detection")
    print(f"   BAM: {args.bam}")
    print(f"   Min depth: {args.min_depth}x")
    
    # Verify if default sites exist in BAM
    samfile = pysam.AlignmentFile(args.bam, "rb")
    bam_refs = samfile.references
    samfile.close()
    
    active_dif_sites = {}
    for site_name, coords in DEFAULT_DIF_SITES.items():
        if coords["chrom"] in bam_refs:
            active_dif_sites[site_name] = coords
        else:
            print(f"   ⚠️  Default site {site_name} ({coords['chrom']}) not in BAM. Attempting BLAST discovery...")
            discovered = discover_region_via_blast(coords["query"], args.reference)
            if discovered:
                print(f"   🎯 Discovered {site_name} at {discovered['chrom']}:{discovered['start']}-{discovered['end']}")
                active_dif_sites[site_name] = discovered
            else:
                print(f"   ❌ Could not locate {site_name} in reference.")
                active_dif_sites[site_name] = None

    # Detect integration
    detection = detect_ctx_integration(Path(args.bam), active_dif_sites, args.min_depth)
    
    # Build report
    report = {
        "sample_id": args.sample or Path(args.bam).stem,
        **detection
    }
    
    # Write output
    with open(args.output, 'w') as f:
        json.dump(report, f, indent=2)
    
    # Print summary
    print("\n📊 Results:")
    print(f"   Status: {detection['ctx_status']}")
    print(f"   Total copy number: {detection['total_copy_number']}")
    
    for site_name, data in detection["integration_sites"].items():
        status_icon = "✅" if data["detected"] else "❌"
        print(f"   {status_icon} {site_name}: {'DETECTED' if data['detected'] else 'NOT DETECTED'}")
        print(f"      Depth: {data['mean_depth']}x, Bridging reads: {data['bridging_reads']}")
    
    if detection["warning"]:
        print(f"\n{detection['warning']}")
    
    print(f"\n📄 Report: {args.output}")


if __name__ == "__main__":
    main()