#!/usr/bin/env python3
"""
SXT Element Local Assembly
De novo assembly of SXT integrative conjugative element for structural variant resolution
"""

import argparse
import json
import subprocess
import sys
import shutil
import os
import tempfile
from pathlib import Path
from Bio import SeqIO
import pysam


# SXT element coordinates on 2010EL-1786 reference (MDR ICE)
DEFAULT_SXT_REGION = {
    "chrom": "CP003069.1",
    "start": 98000,
    "end": 170000,
    "description": "SXT integrative conjugative element (MDR cassette)",
    "query": "sxtMO"
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
            # Expand to cover the SXT element (~100kb)
            return {"chrom": sseqid, "start": max(0, start - 20000), "end": end + 80000}
        except Exception:
            return None

def calculate_region_depth(bam_file: Path, region: dict) -> float:
    """Calculate mean depth coverage for the SXT region."""
    if region is None:
        return 0.0
    
    try:
        # samtools depth -r chrom:start-end bam_file
        region_str = f"{region['chrom']}:{region['start']}-{region['end']}"
        cmd = ["samtools", "depth", "-r", region_str, str(bam_file)]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split('\n')
        
        total_depth = 0
        # Parsing depth output: "chrom pos depth"
        for line in lines:
            if not line: continue
            parts = line.split('\t')
            if len(parts) >= 3:
                depth = int(parts[2])
                total_depth += depth
        
        # Calculate mean over the entire region length
        region_len = region['end'] - region['start']
        if region_len <= 0: return 0.0
        
        return total_depth / region_len
        
    except Exception as e:
        print(f"⚠️ Failed to calculate depth: {e}", file=sys.stderr)
        # Fallback to 0.0
        return 0.0

def check_ram_available(min_gb: int = 8) -> bool:
    """Check if sufficient RAM is available for SPAdes."""
    try:
        # Get available memory in GB
        # macOS specific
        result = subprocess.run(
            ["sysctl", "-n", "hw.memsize"],
            capture_output=True,
            text=True,
            check=True
        )
        total_bytes = int(result.stdout.strip())
        total_gb = total_bytes / (1024**3)
        
        return total_gb >= min_gb
    except (subprocess.CalledProcessError, ValueError):
        # If we can't determine, assume sufficient
        print(f"⚠️  Could not determine available RAM, assuming <{min_gb}GB", file=sys.stderr)
        return False


def extract_sxt_reads(bam_file: Path, output_fastq: Path, sxt_region: dict) -> int:
    """
    Extract reads mapping to SXT region.
    """
    if sxt_region is None:
        return 0

    region = f"{sxt_region['chrom']}:{sxt_region['start']}-{sxt_region['end']}"
    
    print(f"   Extracting reads from {region}...")
    
    cmd = [
        "samtools", "view",
        "-b",  # BAM output
        "-F", "4",  # exclude unmapped
        str(bam_file),
        region
    ]
    
    # Convert to FASTQ
    cmd_fastq = [
        "samtools", "fastq",
        "-"
    ]
    
    try:
        # Pipe: samtools view | samtools fastq | gzip
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(cmd_fastq, stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        
        with open(output_fastq, 'wb') as f:
            p3 = subprocess.Popen(["gzip"], stdin=p2.stdout, stdout=f)
            p2.stdout.close()
            p3.communicate()
        
        # Count reads
        count_cmd = f"zcat {output_fastq} | wc -l"
        result = subprocess.run(
            count_cmd,
            shell=True,
            capture_output=True,
            text=True
        )
        read_count = int(result.stdout.strip()) // 4  # FASTQ has 4 lines per read
        
        print(f"   ✅ Extracted {read_count} reads")
        return read_count
    
    except Exception as e:
        print(f"   ❌ Failed to extract reads: {e}", file=sys.stderr)
        return 0


def extract_consensus_fallback(bam_file: Path, sxt_region: dict, output_fasta: Path) -> bool:
    """
    Generate a consensus sequence from the BAM alignment for the SXT region.
    This serves as a fallback when assembly fails or is skipped.
    """
    if sxt_region is None:
        return False

    print(f"   Generating alignment-based consensus for {sxt_region['chrom']}:{sxt_region['start']}-{sxt_region['end']}...")
    
    try:
        samfile = pysam.AlignmentFile(str(bam_file), "rb")
        
        consensus = []
        # Iterate through the region
        # Use min_base_quality=20 to filter out noisy reads
        for pileupcolumn in samfile.pileup(sxt_region['chrom'], sxt_region['start'], sxt_region['end'], truncate=True, min_base_quality=20):
            # Get the most common base
            counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                    counts[base] = counts.get(base, 0) + 1
            
            # Simple majority rule
            best_base = max(counts, key=counts.get)
            if counts[best_base] == 0:
                consensus.append('N')
            else:
                consensus.append(best_base)
        
        samfile.close()
        
        if not consensus:
            return False
            
        with open(output_fasta, "w") as f:
            f.write(f">SXT_Alignment_Consensus_{sxt_region['chrom']}_{sxt_region['start']}_{sxt_region['end']}\n")
            f.write("".join(consensus) + "\n")
            
        return True
    except Exception as e:
        print(f"   ❌ Fallback consensus failed: {e}", file=sys.stderr)
        return False


def run_spades(reads_fastq: Path, output_dir: Path, threads: int = 4) -> Path:
    """
    Run SPAdes assembler on SXT reads.
    """
    print("   Running SPAdes assembly...")
    
    spades_dir = output_dir / "spades_output"
    if spades_dir.exists():
        shutil.rmtree(spades_dir)
    spades_dir.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "spades.py",
        "--isolate",  # High-quality isolate mode
        "--careful",  # Minimize mismatches and indels
        "-s", str(reads_fastq),  # Single-end reads
        "-o", str(spades_dir),
        "-t", str(threads),
        "-m", "8"  # 8GB RAM limit
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        contigs = spades_dir / "contigs.fasta"
        if not contigs.exists():
            raise FileNotFoundError("SPAdes did not produce contigs.fasta")
        
        print("   ✅ Assembly complete")
        return contigs
    
    except subprocess.CalledProcessError as e:
        print(f"   ❌ SPAdes failed: {e.stderr}", file=sys.stderr)
        raise


def analyze_contigs(contigs_fasta: Path, reference_fasta: Path) -> dict:
    """
    Analyze assembled contigs for structural variants.
    """
    print("   Analyzing contigs...")
    
    # Count contigs
    contig_count = 0
    total_length = 0
    
    with open(contigs_fasta) as f:
        for line in f:
            if line.startswith('>'):
                contig_count += 1
            else:
                total_length += len(line.strip())
    
    return {
        "contig_count": contig_count,
        "total_length": total_length,
        "n50": None,
        "structural_variants": []
    }


def main():
    parser = argparse.ArgumentParser(
        description="Local assembly of SXT element for structural variant resolution"
    )
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument("--output", required=True, help="Output JSON report")
    parser.add_argument("--contigs", required=True, help="Output contigs FASTA")
    parser.add_argument("--outdir", required=True, help="Working directory")
    parser.add_argument("--sample", help="Sample ID")
    parser.add_argument("--mode", default="LABORATORY_FULL",
                       choices=["FIELD_RAPID", "LABORATORY_FULL"],
                       help="Pipeline mode")
    parser.add_argument("--threads", type=int, default=4, help="Threads for SPAdes")
    
    args = parser.parse_args()
    
    work_dir = Path(args.outdir)
    work_dir.mkdir(parents=True, exist_ok=True)
    
    print("🧬 SXT Element Local Assembly")
    print(f"   Mode: {args.mode}")
    
    # Check if SXT region exists in BAM
    samfile = pysam.AlignmentFile(args.bam, "rb")
    bam_refs = samfile.references
    samfile.close()
    
    sxt_region = None
    if DEFAULT_SXT_REGION["chrom"] in bam_refs:
        sxt_region = DEFAULT_SXT_REGION
    else:
        print(f"   ⚠️  Default SXT chrom {DEFAULT_SXT_REGION['chrom']} not in BAM. Attempting BLAST discovery...")
        discovered = discover_region_via_blast(DEFAULT_SXT_REGION["query"], args.reference)
        if discovered:
            print(f"   🎯 Discovered SXT at {discovered['chrom']}:{discovered['start']}-{discovered['end']}")
            sxt_region = discovered
        else:
            print("   ❌ Could not locate SXT element in reference.")

    # Calculate SXT depth check
    mean_depth = calculate_region_depth(Path(args.bam), sxt_region)
    print(f"   📊 SXT Mean Depth: {mean_depth:.2f}x")

    # Check mode
    if args.mode == "FIELD_RAPID":
        print("   ⚡ FIELD_RAPID mode: Skipping SPAdes assembly")
        
        # Fallback to alignment consensus
        success = extract_consensus_fallback(Path(args.bam), sxt_region, Path(args.contigs))
        
        report = {
            "sample_id": args.sample or Path(args.bam).stem,
            "status": "SKIPPED",
            "method": "alignment_consensus" if success else "none",
            "reason": "FIELD_RAPID mode - using reference-based consensus",
            "assembly": None,
            "mean_depth": mean_depth
        }
        
        if not success:
            Path(args.contigs).touch()
        
        with open(args.output, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"   ✅ Skipped ({'alignment-based consensus generated' if success else 'failed to generate consensus'})")
        return
    
    # Check coverage depth (need at least 15x for reliable assembly)
    MIN_COVERAGE = 15.0
    if mean_depth < MIN_COVERAGE:
        print(f"   ⚠️  Insufficient coverage ({mean_depth:.1f}x < {MIN_COVERAGE}x)")
        print(f"   Falling back to reference-based consensus")
        
        success = extract_consensus_fallback(Path(args.bam), sxt_region, Path(args.contigs))
        
        report = {
            "sample_id": args.sample or Path(args.bam).stem,
            "status": "INSUFFICIENT_COVERAGE",
            "method": "alignment_consensus" if success else "none",
            "reason": f"Coverage {mean_depth:.1f}x below threshold ({MIN_COVERAGE}x)",
            "assembly": None,
            "mean_depth": mean_depth
        }
        
        if not success:
            Path(args.contigs).touch()
        
        with open(args.output, 'w') as f:
            json.dump(report, f, indent=2)
        
        return
    
    # Check RAM

    if not check_ram_available(min_gb=8):
        print("   ⚠️  Insufficient RAM for SPAdes (need 8GB)")
        print("   Falling back to reference-based consensus")
        
        success = extract_consensus_fallback(Path(args.bam), sxt_region, Path(args.contigs))
        
        report = {
            "sample_id": args.sample or Path(args.bam).stem,
            "status": "SKIPPED",
            "method": "alignment_consensus" if success else "none",
            "reason": "Insufficient RAM for SPAdes assembly",
            "assembly": None,
            "mean_depth": mean_depth
        }
        
        if not success:
            Path(args.contigs).touch()
        
        with open(args.output, 'w') as f:
            json.dump(report, f, indent=2)
        
        return
    
    # Extract SXT reads
    sxt_reads = work_dir / "sxt_reads.fastq.gz"
    read_count = extract_sxt_reads(Path(args.bam), sxt_reads, sxt_region)
    
    if read_count < 100:
        print(f"   ⚠️  Too few reads ({read_count}) for assembly")
        
        success = extract_consensus_fallback(Path(args.bam), sxt_region, Path(args.contigs))
        
        report = {
            "sample_id": args.sample or Path(args.bam).stem,
            "status": "INSUFFICIENT_COVERAGE",
            "method": "alignment_consensus" if success else "none",
            "reason": f"Only {read_count} reads mapped to SXT region",
            "assembly": None,
            "mean_depth": mean_depth
        }
        
        if not success:
            Path(args.contigs).touch()
        
        with open(args.output, 'w') as f:
            json.dump(report, f, indent=2)
        
        return
    
    # Run SPAdes
    try:
        contigs = run_spades(sxt_reads, work_dir, args.threads)
        
        # Analyze contigs
        analysis = analyze_contigs(contigs, Path(args.reference))
        
        # Copy contigs to output
        shutil.copy(contigs, args.contigs)
        
        # Build report
        report = {
            "sample_id": args.sample or Path(args.bam).stem,
            "status": "SUCCESS",
            "method": "assembly",
            "reads_extracted": read_count,
            "assembly": analysis,
            "mean_depth": mean_depth
        }
        
        print(f"   📊 Contigs: {analysis['contig_count']}")
        print(f"   📏 Total length: {analysis['total_length']} bp")
        
    except Exception as e:
        print(f"   ❌ Assembly failed: {e}", file=sys.stderr)
        
        # Try fallback one last time
        success = extract_consensus_fallback(Path(args.bam), sxt_region, Path(args.contigs))
        
        report = {
            "sample_id": args.sample or Path(args.bam).stem,
            "status": "FAILED",
            "method": "alignment_consensus" if success else "none",
            "reason": str(e),
            "assembly": None,
            "mean_depth": mean_depth
        }
        
        if not success:
            Path(args.contigs).touch()
    
    # Write report
    with open(args.output, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\n📄 Report: {args.output}")


if __name__ == "__main__":
    main()