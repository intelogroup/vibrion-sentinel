#!/usr/bin/env python3
"""
Enhanced Consensus Generation Workflow
Integrates: Platform Detection → Depth Masking → Serogroup Detection

This is the main entry point for consensus generation in the Vibrion pipeline.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path


def run_command(cmd: list, description: str) -> dict:
    """Run a command and return parsed JSON output if applicable."""
    print(f"\n🔧 {description}...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("   ✅ Complete")
        
        # Try to parse JSON output
        try:
            return json.loads(result.stdout)
        except json.JSONDecodeError:
            return {"stdout": result.stdout, "stderr": result.stderr}
    
    except subprocess.CalledProcessError as e:
        print(f"   ❌ Error: {e.stderr}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Enhanced consensus generation with platform detection and QC"
    )
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--fastq", required=True, help="Original FASTQ (for platform detection)")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--mode", default="LABORATORY_FULL", 
                       choices=["FIELD_RAPID", "LABORATORY_FULL"],
                       help="Pipeline mode")
    
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("🧬 VIBRION SENTINEL: Enhanced Consensus Generation")
    print("=" * 60)
    print(f"Sample: {args.sample}")
    print(f"Mode: {args.mode}")
    
    # Step 1: Platform Detection
    platform_report = outdir / "platform_detection.json"
    run_command(
        ["python3", "workflow/scripts/detect_platform.py",
         "--fastq", args.fastq,
         "--output", str(platform_report)],
        "Detecting sequencing platform"
    )
    
    with open(platform_report) as f:
        platform_data = json.load(f)
    
    platform = platform_data["platform"]
    polisher = platform_data["polisher_config"]["tool"]
    
    print(f"\n   Platform: {platform}")
    print(f"   Polisher: {polisher}")
    
    # Step 2: Get mode configuration
    mode_config_result = run_command(
        ["python3", "workflow/scripts/mode_selector.py", args.mode],
        f"Loading {args.mode} configuration"
    )
    
    min_depth = mode_config_result["min_depth"]
    print(f"   Min depth threshold: {min_depth}x")
    
    # Step 3: Generate consensus with depth masking
    print(f"\n{'=' * 60}")
    print("📖 Generating consensus genome...")
    print(f"{'=' * 60}")
    
    consensus_result = run_command(
        ["python3", "workflow/scripts/generate_consensus_genome.py",
         "--bam", args.bam,
         "--ref", args.ref,
         "--outdir", str(outdir),
         "--sample", args.sample,
         "--min-depth", str(min_depth),
         "--min-qual", "20"],
        "Running consensus caller"
    )
    
    # Step 4: Serogroup Detection (O1 vs O139 vs Phage Scar)
    serogroup_report = outdir / "serogroup_detection.json"
    run_command(
        ["python3", "workflow/scripts/detect_serogroup.py",
         "--bam", args.bam,
         "--reference", args.ref,
         "--output", str(serogroup_report)],
        "Detecting serogroup (O1/O139/Phage Scar)"
    )
    
    with open(serogroup_report) as f:
        serogroup_data = json.load(f)
    
    print(f"\n   Serogroup: {serogroup_data['serogroup']}")
    print(f"   Confidence: {serogroup_data['confidence']}")
    print(f"   Reason: {serogroup_data['reason']}")
    
    # Step 5: Generate final integrated report
    final_report = {
        "sample_id": args.sample,
        "mode": args.mode,
        "platform": platform_data,
        "serogroup": serogroup_data,
        "consensus_stats": consensus_result.get("statistics", {}),
        "output_files": {
            "consensus_fasta": str(outdir / f"{args.sample}_consensus.fasta"),
            "platform_report": str(platform_report),
            "serogroup_report": str(serogroup_report)
        },
        "warnings": []
    }
    
    # Add warnings based on results
    if serogroup_data["serogroup"] == "PHAGE_SCAR_POSSIBLE":
        final_report["warnings"].append(
            "⚠️ Possible phage predation damage detected. rfb region has low coverage without bridging reads."
        )
    
    if serogroup_data["serogroup"] == "O139_CANDIDATE":
        final_report["warnings"].append(
            "🚨 Potential O139 serogroup detected. Requires wbf marker confirmation."
        )
    
    # Check for polyclonal infection
    het_sites = consensus_result.get("statistics", {}).get("heterogeneous_sites_count", 0)
    if het_sites > 50:
        final_report["warnings"].append(
            f"⚠️ High heterozygosity detected ({het_sites} sites). Possible polyclonal infection."
        )
    
    # Save final report
    final_report_path = outdir / "consensus_generation_summary.json"
    with open(final_report_path, 'w') as f:
        json.dump(final_report, f, indent=2)
    
    print(f"\n{'=' * 60}")
    print("✅ CONSENSUS GENERATION COMPLETE")
    print(f"{'=' * 60}")
    print(f"📁 Output directory: {outdir}")
    print(f"📄 Final report: {final_report_path}")
    
    if final_report["warnings"]:
        print("\n⚠️  WARNINGS:")
        for warning in final_report["warnings"]:
            print(f"   {warning}")
    
    print("\n🎯 Next steps: Variant calling → SnpEff annotation → Evo2 analysis")


if __name__ == "__main__":
    main()
