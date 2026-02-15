import subprocess
import os
from pathlib import Path

def forensic_serotype_audit(sample_path):
    print(f"🕵️ MULTI-LINEAGE DIVERGENCE AUDIT: {sample_path}")
    print("=" * 80)
    
    # 1. Signature Definitions
    # Haitian Ogawa (Parent/Resilient): Functional wbeT sequence
    sig_ogawa = "TTTTAAGTGA"
    # Haitian Inaba (Adaptive/2015): 1-bp deletion at nt-305
    sig_inaba = "TTTAAGTGA"
    # SXT Haiti Signature (Clonal Baseline)
    sig_sxt = "GATTGGCGGTTTT"
    
    def count_hits(sig):
        cmd = f"gzcat {sample_path} | head -n 500000 | grep -c {sig}"
        try:
            res = subprocess.check_output(cmd, shell=True).decode().strip()
            return int(res)
        except:
            return 0

    print("🚀 Scanning raw reads for lineage markers (500k slice)...")
    hits_ogawa = count_hits(sig_ogawa)
    hits_inaba = count_hits(sig_inaba)
    hits_sxt = count_hits(sig_sxt)
    
    # 2. Forensic Interpretation
    is_haitian_clone = hits_sxt > 0
    
    print(f"\n📊 RESULTS:")
    print(f"   Haitian SXT (Baseline): {hits_sxt} hits")
    print(f"   Ogawa (Parent):         {hits_ogawa} hits")
    print(f"   Inaba (Adaptive):       {hits_inaba} hits")
    
    print(f"\n🧠 PIPELINE INTELLIGENCE:")
    if not is_haitian_clone:
        print("   🔴 VERDICT: NON-HAITIAN LINEAGE (Potential New Introduction)")
    else:
        if hits_ogawa > 0 and hits_inaba > 0:
            print("   🟡 VERDICT: MULTI-LINEAGE CO-PERSISTENCE (Mixed Population)")
            print("      Note: Both Ogawa and Inaba branches detected. High risk of immunological escape.")
        elif hits_ogawa > hits_inaba:
            print("   🟢 VERDICT: HAITIAN OGAWA (Clonal Continuation)")
            print("      Note: This lineage fueled the 2022 reemergence.")
        elif hits_inaba > hits_ogawa:
            print("   🔵 VERDICT: HAITIAN INABA (2015 Adaptive Branch)")
            print("      Note: Dominated 2015-2019; currently co-persisting.")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        forensic_serotype_audit(sys.argv[1])
    else:
        # Test on the 2021 environmental sample
        forensic_serotype_audit("data/raw_reads/SRR22265446_1.fastq.gz")
