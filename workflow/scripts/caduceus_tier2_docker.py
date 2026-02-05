#!/usr/bin/env python3
"""
Tier 2: Caduceus Gene-Level Inspection (Docker Version)

This script runs REAL Caduceus model with mamba-ssm in a Docker container.
Designed to be called from the host system via Docker.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

import torch
from Bio import SeqIO

# Import Caduceus dependencies
try:
    from transformers import AutoTokenizer, AutoModel
    import mamba_ssm
    print("✅ mamba-ssm available", file=sys.stderr)
except ImportError as e:
    print(f"❌ Error: {e}", file=sys.stderr)
    print("This script requires mamba-ssm. Run in Docker container.", file=sys.stderr)
    sys.exit(1)


def load_fasta(fasta_path: Path) -> Dict[str, str]:
    """Load sequences from FASTA file."""
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def tier2_caduceus_inspect(
    sample_seq: str,
    loci_seqs: Dict[str, str],
    model_id: str,
    cache_dir: str,
    device: torch.device,
    threshold: float,
    anomaly_windows: List[Dict],
    target_loci: List[str],
    overlap_threshold: float
) -> Dict:
    """
    Tier 2: Caduceus Gene-Level Inspection with REAL mamba-ssm
    """
    print("\n🧬 Tier 2: Caduceus Gene-Level Inspection (Docker)")
    print(f"   Model: {model_id}")
    print(f"   Device: {device}")
    
    try:
        print("   → Loading Caduceus tokenizer...")
        tokenizer = AutoTokenizer.from_pretrained(
            model_id,
            cache_dir=cache_dir,
            trust_remote_code=True
        )
        
        print("   → Loading Caduceus model...")
        model = AutoModel.from_pretrained(
            model_id,
            cache_dir=cache_dir,
            trust_remote_code=True,
            torch_dtype=torch.float32
        )
        model = model.to(device)
        model.eval()
        
        print("   ✅ Caduceus model loaded with mamba-ssm")
        
        # Analyze target loci
        locus_scores = {}
        high_risk_loci = []
        
        with torch.no_grad():
            for locus_name in target_loci:
                if locus_name not in loci_seqs:
                    print(f"   ⚠️  Locus {locus_name} not found")
                    continue
                
                locus_seq = loci_seqs[locus_name]
                print(f"   → Analyzing {locus_name} ({len(locus_seq)} bp)...")
                
                # Tokenize and run inference
                inputs = tokenizer(locus_seq, return_tensors="pt", truncation=True, max_length=1024)
                inputs = {k: v.to(device) for k, v in inputs.items()}
                
                outputs = model(**inputs)
                
                # Compute risk score from model outputs
                if hasattr(outputs, 'last_hidden_state'):
                    risk_score = torch.norm(outputs.last_hidden_state, dim=-1).mean().item()
                else:
                    risk_score = 0.5
                
                locus_scores[locus_name] = round(risk_score, 4)
                
                if risk_score > threshold:
                    high_risk_loci.append(locus_name)
                    print(f"      ⚠️  HIGH RISK: {risk_score:.4f}")
                else:
                    print(f"      ✅ Safe: {risk_score:.4f}")
        
        escalate = len(high_risk_loci) > 0
        
        result = {
            "tier": 2,
            "method": "Caduceus (Docker + mamba-ssm)",
            "model": model_id,
            "locus_scores": locus_scores,
            "high_risk_loci": high_risk_loci,
            "threshold": threshold,
            "escalate_to_cloud": escalate,
            "decision": f"ESCALATE - {len(high_risk_loci)} high-risk loci" if escalate else "LOW RISK"
        }
        
        print("\n   📊 Results:")
        print(f"      Locus Scores: {locus_scores}")
        print(f"      High Risk: {high_risk_loci}")
        print(f"      Decision: {result['decision']}")
        
        return result
        
    except Exception as e:
        print(f"   ❌ Error: {e}", file=sys.stderr)
        return {
            "tier": 2,
            "method": "Caduceus",
            "model": model_id,
            "error": str(e),
            "escalate_to_cloud": True,
            "decision": "ERROR - Escalating to be safe"
        }


def main():
    parser = argparse.ArgumentParser(description="Tier 2: Caduceus Docker Analysis")
    parser.add_argument("--tier1-json", required=True, help="Tier 1 output JSON")
    parser.add_argument("--sample-fasta", required=True, help="Sample consensus FASTA")
    parser.add_argument("--loci-fasta", required=True, help="Surveillance loci FASTA")
    parser.add_argument("--model", required=True, help="Caduceus model ID")
    parser.add_argument("--loci", nargs="+", default=["ctxB", "wbeT", "tcpA"], help="Target loci")
    parser.add_argument("--threshold", type=float, default=0.5, help="Risk threshold")
    parser.add_argument("--overlap-threshold", type=float, default=0.5, help="Overlap threshold")
    parser.add_argument("--device", default="cuda", help="Device: cuda or cpu")
    parser.add_argument("--cache-dir", default="/workspace/data/models/triage", help="Model cache")
    parser.add_argument("--output", required=True, help="Output JSON path")
    
    args = parser.parse_args()
    
    print("🐳 Caduceus Tier 2 Analysis (Docker + mamba-ssm)")
    print("=" * 60)
    
    # Load Tier 1 results
    with open(args.tier1_json) as f:
        tier1_result = json.load(f)
    
    # Load sequences
    sample_seqs = load_fasta(Path(args.sample_fasta))
    loci_seqs = load_fasta(Path(args.loci_fasta))
    
    sample_seq = list(sample_seqs.values())[0] if sample_seqs else ""
    
    # Get device
    if args.device == "cuda" and torch.cuda.is_available():
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")
    
    # Run Tier 2
    result = tier2_caduceus_inspect(
        sample_seq=sample_seq,
        loci_seqs=loci_seqs,
        model_id=args.model,
        cache_dir=args.cache_dir,
        device=device,
        threshold=args.threshold,
        anomaly_windows=tier1_result.get("anomaly_windows", []),
        target_loci=args.loci,
        overlap_threshold=args.overlap_threshold
    )
    
    # Write output
    with open(args.output, 'w') as f:
        json.dump(result, f, indent=2)
    
    print(f"\n📄 Results written to: {args.output}")
    return 0


if __name__ == "__main__":
    exit(main())
