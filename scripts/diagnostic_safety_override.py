import sys

def diagnostic_override(fastq_path):
    print(f"🩺 DIAGNOSTIC SAFETY OVERRIDE ANALYSIS: {fastq_path}")
    print("=" * 80)
    
    # 1. Count Phage vs Vibrio (Using our simulated signatures)
    phage_count = 0
    vibrio_count = 0
    
    with open(fastq_path, 'r') as f:
        for line in f:
            if "ACGAAGTGGAACGC" in line:
                phage_count += 1
            if "GCAAAGGCGATTCG" in line:
                vibrio_count += 1
                
    ratio = phage_count / max(vibrio_count, 1)
    print(f"🔍 Phage Count: {phage_count}")
    print(f"🔍 Vibrio Count: {vibrio_count}")
    print(f"📊 Phage/Vibrio Ratio: {ratio:.2f}")
    
    # 2. Decision Logic (Threshold 0.05)
    if ratio > 0.05:
        print("\n🚨 CRITICAL ALERT: HIGH PHAGE PREDATION DETECTED")
        print("   Status: CULTURE RISK HIGH (Lytic phage activity likely killing Vibrio on plates)")
        print("   Recommendation: OVERRIDE CULTURE. PROCEED WITH PCR-FIRST (ctxAB, rfb).")
        print("   Reasoning: ICP2-like phage abundance suggests traditional culture will yield false negatives.")
    else:
        print("\n✅ DIAGNOSTIC CONFIDENCE: HIGH")
        print("   Status: Phage levels low. Culture results reliable.")

if __name__ == "__main__":
    diagnostic_override("data/tests/phage_test/simulated_phage.fastq")
