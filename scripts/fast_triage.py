#!/usr/bin/env python3
"""
Vibrion Sentinel - Fast Triage System
=====================================
30-second triage for cholera samples using k-mer matching.

Usage:
    python fast_triage.py --input sample.fastq.gz --output report.json
    python fast_triage.py --input /data/input --output /data/output --batch

Features:
    - K-mer based lineage matching (Sourmash)
    - AMR gene k-mer scanning
    - Risk scoring
    - Conditional escalation to full pipeline
"""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from datetime import datetime

# Try to import optional dependencies
try:
    from rich.console import Console
    from rich.table import Table
    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False


class FastTriage:
    """High-speed k-mer based triage for cholera surveillance."""
    
    # Confidence thresholds
    HIGH_CONFIDENCE_THRESHOLD = 0.85  # Skip full pipeline
    ANOMALY_THRESHOLD = 0.50  # Definitely needs full pipeline
    
    # Critical AMR markers that auto-escalate risk
    CRITICAL_AMR = {"blaNDM-1", "erm(B)", "mph(E)", "mph(A)", "gyrA", "qnrS", "qnrB"}
    
    def __init__(self, references_dir=None, verbose=True):
        self.verbose = verbose
        self.references_dir = Path(references_dir or os.environ.get(
            "VIBRION_REFERENCES", 
            Path(__file__).parent.parent / "data" / "global_references"
        ))
        self.console = Console() if RICH_AVAILABLE else None
        
    def log(self, msg, style=""):
        """Print status message."""
        if self.verbose:
            if self.console:
                self.console.print(f"[{style}]{msg}[/{style}]" if style else msg)
            else:
                print(msg)
    
    def run_sourmash_sketch(self, input_file, output_sig):
        """Generate k-mer signature for input file."""
        cmd = [
            "sourmash", "sketch", "dna",
            "-p", "k=31,scaled=1000",
            str(input_file),
            "-o", str(output_sig),
            "--name", input_file.stem
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result.returncode == 0
    
    def run_sourmash_search(self, sample_sig, reference_sigs):
        """Compare sample against reference signatures."""
        if not reference_sigs:
            return []
        
        # Combine reference sigs into temp file list
        sig_files = " ".join(str(s) for s in reference_sigs)
        cmd = f"sourmash search -q {sample_sig} {sig_files} --threshold 0.1 -o /dev/stdout"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        matches = []
        if result.returncode == 0:
            lines = result.stdout.strip().split('\n')
            if len(lines) > 1:  # Has header + data
                for line in lines[1:]:
                    parts = line.split(',')
                    if len(parts) >= 4:
                        matches.append({
                            "similarity": float(parts[0]),
                            "name": parts[3].strip('"')
                        })
        return sorted(matches, key=lambda x: x["similarity"], reverse=True)
    
    def scan_amr_genes(self, input_file, output_dir):
        """Scan for AMR genes using KMA alignment against CARD database."""
        detected = []
        
        # Path to CARD KMA index
        amr_db = Path(__file__).parent.parent / "data" / "amr_db" / "card_index"
        
        if not amr_db.with_suffix(".name").exists():
            self.log(f"⚠️  KMA index not found at {amr_db}", "yellow")
            return []
        
        # Run KMA alignment
        output_prefix = output_dir / "kma_amr"
        cmd = [
            "kma",
            "-i", str(input_file),
            "-o", str(output_prefix),
            "-t_db", str(amr_db),
            "-t", "4",           # threads
            "-1t1",              # one template per query
            "-cge",              # CGE output format
            "-nf",               # no fragments
            "-mem_mode",         # memory efficient mode
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Parse KMA results (.res file)
        res_file = Path(f"{output_prefix}.res")
        if res_file.exists():
            import re
            with open(res_file) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 1:
                        # Extract gene name from template name
                        # Format: gb|XXX|+|position|ARO:XXX|GeneName [Species]
                        template = parts[0]
                        # Try to extract gene name
                        gene_match = re.search(r'\|([A-Za-z]+[\(\)A-Za-z0-9\-]+)\s*\[', template)
                        if gene_match:
                            gene_name = gene_match.group(1)
                            if gene_name not in detected:
                                detected.append(gene_name)
                        else:
                            # Fallback: try different pattern
                            match2 = re.search(r'ARO:\d+\|([^\s\[]+)', template)
                            if match2:
                                gene_name = match2.group(1)
                                if gene_name not in detected:
                                    detected.append(gene_name)
        
        return detected

    
    def calculate_risk_score(self, lineage_match, amr_detected):
        """Calculate overall risk score (0-100)."""
        score = 0
        
        # Base score from lineage match
        if lineage_match:
            best_sim = lineage_match[0]["similarity"]
            # Known toxigenic lineages get higher base score
            toxigenic_lineages = {"haiti-2010", "haiti-2022", "yemen", "malawi", "bangladesh", "kenya"}
            if any(tox in lineage_match[0]["name"].lower() for tox in toxigenic_lineages):
                score += 60 * best_sim
            else:
                score += 20 * best_sim
        
        # AMR contribution
        for gene in amr_detected:
            if gene in self.CRITICAL_AMR:
                score += 15  # Critical resistance
            else:
                score += 5   # Standard resistance
        
        return min(100, int(score))
    
    def determine_escalation(self, confidence, amr_detected, risk_score):
        """Determine if full pipeline is needed."""
        # Auto-escalate if:
        # 1. Low confidence match
        # 2. Critical AMR detected
        # 3. High risk score
        
        if confidence < self.ANOMALY_THRESHOLD:
            return True, "LOW_CONFIDENCE"
        
        if any(gene in self.CRITICAL_AMR for gene in amr_detected):
            return True, "CRITICAL_AMR"
        
        if risk_score >= 70:
            return True, "HIGH_RISK"
        
        if confidence < self.HIGH_CONFIDENCE_THRESHOLD:
            return True, "MODERATE_CONFIDENCE"
        
        return False, "HIGH_CONFIDENCE"
    
    def triage(self, input_file, output_dir=None):
        """Run fast triage on a single sample."""
        input_path = Path(input_file)
        output_dir = Path(output_dir or input_path.parent)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        start_time = time.time()
        sample_name = input_path.stem.replace("_1", "").replace("_2", "").replace(".fastq", "")
        
        self.log(f"\n🧬 Vibrion Fast Triage: {sample_name}", "bold cyan")
        self.log("=" * 50, "dim")
        
        # Step 1: Generate k-mer signature
        self.log("📊 Generating k-mer signature...", "yellow")
        sig_file = output_dir / f"{sample_name}.sig"
        if not self.run_sourmash_sketch(input_path, sig_file):
            return {"error": "Failed to generate k-mer signature"}
        
        # Step 2: Compare to reference library
        self.log("🔍 Matching against reference library...", "yellow")
        reference_sigs = list(self.references_dir.glob("*.sig"))
        matches = self.run_sourmash_search(sig_file, reference_sigs)
        
        confidence = matches[0]["similarity"] if matches else 0.0
        best_match = matches[0]["name"] if matches else "Unknown"
        
        # Step 3: Scan for AMR genes using KMA
        self.log("💊 Scanning for AMR markers (KMA)...", "yellow")
        amr_detected = self.scan_amr_genes(input_path, output_dir)
        
        # Step 4: Calculate risk score
        risk_score = self.calculate_risk_score(matches, amr_detected)
        
        # Step 5: Determine escalation
        needs_full_pipeline, reason = self.determine_escalation(
            confidence, amr_detected, risk_score
        )
        
        elapsed = time.time() - start_time
        
        # Build result
        result = {
            "sample": sample_name,
            "timestamp": datetime.now().isoformat(),
            "elapsed_seconds": round(elapsed, 2),
            "lineage": {
                "best_match": best_match,
                "confidence": round(confidence, 4),
                "all_matches": matches[:5]
            },
            "amr": {
                "detected": amr_detected,
                "critical": [g for g in amr_detected if g in self.CRITICAL_AMR]
            },
            "risk_score": risk_score,
            "escalation": {
                "needs_full_pipeline": needs_full_pipeline,
                "reason": reason
            },
            "verdict": self._generate_verdict(confidence, amr_detected, risk_score)
        }
        
        # Save JSON report
        json_file = output_dir / f"{sample_name}_triage.json"
        with open(json_file, 'w') as f:
            json.dump(result, f, indent=2)
        
        # Save Markdown report
        md_file = output_dir / f"{sample_name}_triage.md"
        self._write_markdown_report(result, md_file)
        
        # Print summary
        self._print_summary(result)
        
        return result
    
    def _generate_verdict(self, confidence, amr_detected, risk_score):
        """Generate human-readable verdict."""
        if risk_score >= 80:
            return "🔴 HIGH RISK - Immediate action required"
        elif risk_score >= 50:
            return "🟠 MODERATE RISK - Enhanced surveillance"
        elif risk_score >= 25:
            return "🟡 LOW RISK - Standard monitoring"
        else:
            return "🟢 MINIMAL RISK - Routine surveillance"
    
    def _print_summary(self, result):
        """Print formatted summary."""
        self.log("\n" + "=" * 50, "dim")
        self.log(f"⏱️  Completed in {result['elapsed_seconds']}s", "green")
        self.log(f"🧬 Lineage: {result['lineage']['best_match']} ({result['lineage']['confidence']:.1%})", "cyan")
        
        if result['amr']['detected']:
            amr_str = ", ".join(result['amr']['detected'])
            style = "bold red" if result['amr']['critical'] else "yellow"
            self.log(f"💊 AMR: {amr_str}", style)
        else:
            self.log("💊 AMR: None detected", "green")
        
        self.log(f"📈 Risk Score: {result['risk_score']}/100", "bold")
        self.log(f"\n{result['verdict']}", "bold")
        
        if result['escalation']['needs_full_pipeline']:
            self.log(f"\n⚠️  Full pipeline recommended: {result['escalation']['reason']}", "bold yellow")
    
    def _write_markdown_report(self, result, output_file):
        """Write simplified Markdown report."""
        md = f"""# Vibrion Fast Triage Report

**Sample:** {result['sample']}  
**Date:** {result['timestamp']}  
**Processing Time:** {result['elapsed_seconds']}s

---

## Lineage Match
- **Best Match:** {result['lineage']['best_match']}
- **Confidence:** {result['lineage']['confidence']:.1%}

## AMR Detection
- **Detected:** {', '.join(result['amr']['detected']) or 'None'}
- **Critical:** {', '.join(result['amr']['critical']) or 'None'}

## Risk Assessment
- **Score:** {result['risk_score']}/100
- **Verdict:** {result['verdict']}

## Escalation
- **Full Pipeline Needed:** {'Yes' if result['escalation']['needs_full_pipeline'] else 'No'}
- **Reason:** {result['escalation']['reason']}

---
*Generated by Vibrion Sentinel Fast Triage*
"""
        with open(output_file, 'w') as f:
            f.write(md)


def main():
    parser = argparse.ArgumentParser(
        description="Vibrion Sentinel - Fast Triage (30-second cholera sample analysis)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Single sample:
    python fast_triage.py --input sample.fastq.gz --output ./results

  Batch processing:
    python fast_triage.py --input /data/samples --output /data/results --batch

  With custom references:
    python fast_triage.py --input sample.fastq.gz --references /path/to/refs
        """
    )
    parser.add_argument("--input", "-i", required=True, help="Input FASTQ file or directory")
    parser.add_argument("--output", "-o", help="Output directory (default: same as input)")
    parser.add_argument("--references", "-r", help="Reference signatures directory")
    parser.add_argument("--batch", "-b", action="store_true", help="Process all files in input directory")
    parser.add_argument("--quiet", "-q", action="store_true", help="Suppress progress output")
    parser.add_argument("--json", action="store_true", help="Output only JSON to stdout")
    
    args = parser.parse_args()
    
    triage = FastTriage(references_dir=args.references, verbose=not args.quiet and not args.json)
    
    input_path = Path(args.input)
    
    if args.batch and input_path.is_dir():
        # Batch mode
        results = []
        for fastq in input_path.glob("*.fastq*"):
            result = triage.triage(fastq, args.output)
            results.append(result)
        
        if args.json:
            print(json.dumps(results, indent=2))
    else:
        # Single file mode
        result = triage.triage(input_path, args.output)
        
        if args.json:
            print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
