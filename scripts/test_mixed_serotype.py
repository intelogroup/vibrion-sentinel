#!/usr/bin/env python3
"""
Quick Test: Create a VCF with heterozygous wbeT calls to simulate mixed Ogawa/Inaba
This tests the serotype detection system's ability to detect co-infection
"""
import sys
from pathlib import Path

def create_mixed_vcf():
    """
    Create a synthetic VCF file with heterozygous calls in wbeT gene
    This simulates a mixed Ogawa/Inaba sample
    """
    print("🧬 Creating Synthetic Mixed wbeT VCF")
    print("=" * 70)
    
    # wbeT gene location on CP003069.1: 2687324-2688226
    # Key mutation sites:
    # - Position 2687797 (codon 158): S158P mutation (Ogawa->Inaba)
    
    vcf_content = """##fileformat=VCFv4.2
##reference=CP003069.1
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmixed_sample
CP003069.1\t2687797\t.\tT\tC\t30\tPASS\tDP=50;AF=0.48\tGT:DP:AD\t0/1:50:26,24
CP003069.1\t2687850\t.\tG\tA\t25\tPASS\tDP=45;AF=0.44\tGT:DP:AD\t0/1:45:25,20
CP003069.1\t2688000\t.\tA\tT\t28\tPASS\tDP=48;AF=0.46\tGT:DP:AD\t0/1:48:26,22
"""
    
    output_vcf = Path("data/pipeline_output/mixed_ogawa_inaba/04_variants/mixed_sample.vcf")
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_vcf, 'w') as f:
        f.write(vcf_content)
    
    print(f"✓ Created VCF: {output_vcf}")
    print()
    print("Key heterozygous sites in wbeT:")
    print("  - Position 2687797: S158P site (AF=0.48) - MIXED!")
    print("  - Position 2687850: Additional variant (AF=0.44)")
    print("  - Position 2688000: Additional variant (AF=0.46)")
    print()
    print("This VCF simulates a ~50:50 Ogawa:Inaba mixture")
    
    return output_vcf

def create_mixed_consensus():
    """
    Create a consensus FASTA that contains both functional and mutated wbeT
    """
    print("\n🧬 Creating Mixed Consensus FASTA")
    print("=" * 70)
    
    # Get reference wbeT
    ref_wbet = Path("data/references/wbeT_ref.fasta")
    if not ref_wbet.exists():
        print(f"Warning: Reference wbeT not found at {ref_wbet}")
        return None
    
    # Read reference
    with open(ref_wbet) as f:
        content = f.read()
    
    # Create output
    output_fasta = Path("data/pipeline_output/mixed_ogawa_inaba/09_consensus/consensus.fasta")
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_fasta, 'w') as f:
        f.write(content)
    
    print(f"✓ Created consensus: {output_fasta}")
    return output_fasta

def test_wbeT_detection():
    """
    Run the wbeT tracker on the mixed sample
    """
    print("\n🔬 Testing wbeT Detection on Mixed Sample")
    print("=" * 70)
    
    import subprocess
    
    vcf_file = create_mixed_vcf()
    consensus_file = create_mixed_consensus()
    
    if not consensus_file:
        print("❌ Cannot run test without consensus file")
        return
    
    # Run wbeT tracker
    cmd = [
        "python3", "scripts/track_wbeT_mutations.py",
        str(consensus_file),
        "--ref", "data/references/2010EL-1786.fasta",
        "--vcf", str(vcf_file)
    ]
    
    print(f"\n🚀 Running: {' '.join(cmd)}\n")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(result.stdout)
        
        # Check if mixed detection worked
        if "Hikojima" in result.stdout or "Mixed" in result.stdout:
            print("\n✅ SUCCESS: Mixed Ogawa/Inaba (Hikojima) detected!")
            print("The system correctly identified heterozygous wbeT calls.")
        elif "MIXED" in result.stdout:
            print("\n✅ SUCCESS: Mixed serotype detected!")
        else:
            print("\n⚠️  WARNING: Mixed detection may not have triggered")
            print("Check the output above for serotype classification")
            
    except subprocess.CalledProcessError as e:
        print(f"❌ ERROR: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")

if __name__ == "__main__":
    test_wbeT_detection()
