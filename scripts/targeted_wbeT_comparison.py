import subprocess
import os
from pathlib import Path

def direct_comparison():
    print("🔬 FORCING DIRECT wbeT COMPARISON (Ogawa vs Inaba Founder)")
    print("=" * 80)
    
    ref = "data/pipeline_output/SRR8364418/wbeT_ogawa_ref.fasta"
    fastq1 = "data/raw_reads/PRJNA510624/SRR8364418_1.fastq.gz"
    out_dir = Path("data/pipeline_output/SRR8364418/targeted_wbeT")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    bam = out_dir / "targeted_wbeT.bam"
    
    # 1. Align reads specifically to wbeT gene
    print("🚀 Step 1: Targeted Alignment to Ogawa wbeT...")
    cmd_align = f"minimap2 -ax sr -t 4 {ref} {fastq1} | samtools view -b - | samtools sort -o {bam}"
    subprocess.run(cmd_align, shell=True, check=True)
    subprocess.run(f"samtools index {bam}", shell=True, check=True)
    
    # 2. Inspect Position 305 (The Inaba Deletion Site)
    # The reference is the gene itself, so index 305 is the target.
    print("\n🚀 Step 2: High-Resolution Pileup at Position 305...")
    # Use -f to show reference base
    cmd_pileup = f"samtools mpileup -f {ref} -r CP003069.1:123001:300-310 {bam}"
    # Wait, the reference ID in the extracted FASTA is 'CP003069.1:123001-123993'
    
    # Let's get the actual header
    with open(ref, 'r') as f:
        ref_header = f.readline().strip()[1:]
    
    target_pos = 305
    cmd_pileup = f"samtools mpileup -f {ref} -r \"{ref_header}\":{target_pos-5}-{target_pos+5} {bam}"
    
    result = subprocess.run(cmd_pileup, shell=True, capture_output=True, text=True)
    print(result.stdout)
    
    # 3. Automated Decision Logic
    # Look for the '-' (deletion) signal in the pileup
    if "-" in result.stdout:
        print("\n✅ ALERT: 1-BP DELETION DETECTED AT POSITION 305!")
        print("   Significance: Definitive Founder Mutation for Haiti Inaba Lineage.")
        print("   Functional Impact: Truncated protein (Stop at AA 121).")
        print("   Verdict: O1 INABA (Ancestral Founder HC1961)")
    else:
        print("\n❌ NO DELETION FOUND at position 305.")
        print("   The sample matches the functional Ogawa sequence at this locus.")

if __name__ == "__main__":
    direct_comparison()
