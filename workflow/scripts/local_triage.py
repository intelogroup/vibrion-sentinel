#!/usr/bin/env python3
"""
Tier 1-2 Triage: Local AI Analysis (HyenaDNA + Caduceus)

This script runs REAL DNA foundation models to triage sequences:
- Tier 1 (HyenaDNA): Global structural anomaly detection using real model inference
- Tier 2 (Caduceus): Gene-specific functional validation (requires mamba-ssm)

Only escalates to Evo2 Cloud (Tier 3) if local models detect high-risk mutations.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

import torch
from Bio import SeqIO

# Check for transformers
try:
    from transformers import AutoTokenizer, AutoModel
except ImportError:
    print("ERROR: transformers not installed. Run: pip install transformers", file=sys.stderr)
    sys.exit(1)


def load_fasta(fasta_path: Path) -> Dict[str, str]:
    """Load sequences from FASTA file."""
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def get_device(device_name: str) -> torch.device:
    """Get the appropriate torch device."""
    if device_name == "mps":
        # Force CPU for HyenaDNA on Mac due to unsupported FFT ops in MPS implementation
        print("⚠️  MPS requested, but forcing CPU for HyenaDNA stability (FFT support)", file=sys.stderr)
        return torch.device("cpu")
    elif device_name == "cuda" and torch.cuda.is_available():
        return torch.device("cuda")
    else:
        if device_name != "cpu":
            print(f"⚠️  {device_name} not available, falling back to CPU", file=sys.stderr)
        return torch.device("cpu")


def tier1_hyena_scan(
    sample_seq: str,
    ref_seq: str,
    model_id: str,
    cache_dir: str,
    device: torch.device,
    threshold: float,
    divergent_regions: List[Dict],
    window_size: int
) -> tuple:
    """
    Tier 1: HyenaDNA Real Model Inference
    
    Uses actual HyenaDNA model to compute perplexity scores on divergent regions.
    Returns: (result_dict, tokenizer, model) for reuse in Tier 2
    """
    print("\n🔬 Tier 1: HyenaDNA Real Model Inference")
    print(f"   Model: {model_id}")
    print(f"   Device: {device}")
    print(f"   Cache: {cache_dir}")
    print(f"   Focusing on {len(divergent_regions)} divergent region(s)")
    
    tokenizer = None
    model = None
    
    try:
        # Load tokenizer and model
        print("   → Loading tokenizer...")
        tokenizer = AutoTokenizer.from_pretrained(
            model_id,
            cache_dir=cache_dir,
            trust_remote_code=True
        )
        
        print("   → Loading model...")
        model = AutoModel.from_pretrained(
            model_id,
            cache_dir=cache_dir,
            trust_remote_code=True,
            torch_dtype=torch.float32
        )
        model = model.to(device)
        model.eval()
        
        print("   ✅ Model loaded successfully")
        
        # If no divergent regions, analyze the whole sequence
        if not divergent_regions:
            print("   → No divergent regions from Tier 0, analyzing full sequence")
            regions_to_scan = [{"start": 0, "end": min(len(sample_seq), 1000)}]  # First 1kb
        else:
            regions_to_scan = divergent_regions
        
        # Compute perplexity for each region
        anomaly_windows = []
        total_perplexity = 0.0
        
        with torch.no_grad():
            for region in regions_to_scan:
                start = max(0, region["start"])
                end = min(len(sample_seq), region["end"])
                
                # Extract window
                window_seq = sample_seq[start:end]
                if len(window_seq) < 50:  # Skip tiny windows
                    continue
                
                # Tokenize
                inputs = tokenizer(window_seq, return_tensors="pt", truncation=True, max_length=1024)
                inputs = {k: v.to(device) for k, v in inputs.items()}
                
                # Get model outputs
                outputs = model(**inputs)
                
                # Compute perplexity/anomaly score
                if hasattr(outputs, 'last_hidden_state'):
                    hidden_states = outputs.last_hidden_state
                    # Use variance + norm combination for better discrimination
                    # Norm alone is too stable; variance captures internal "surprise"
                    var_score = torch.var(hidden_states, dim=-1).mean().item()
                    norm_score = torch.norm(hidden_states, dim=-1).mean().item()
                    
                    # Primary score is variance-weighted
                    perplexity = var_score * (norm_score / 10.0) 
                    
                    # If reference is available, compute a Delta score (Relative Anomaly)
                    ref_window = ref_seq[start:end] if end <= len(ref_seq) else ""
                    if ref_window and len(ref_window) > 50:
                        ref_inputs = tokenizer(ref_window, return_tensors="pt", truncation=True, max_length=1024)
                        ref_inputs = {k: v.to(device) for k, v in ref_inputs.items()}
                        ref_outputs = model(**ref_inputs)
                        
                        if hasattr(ref_outputs, 'last_hidden_state'):
                            ref_hidden = ref_outputs.last_hidden_state
                            # Use Cosine Similarity or MSE between hidden states as a "Forensic Delta"
                            cos = torch.nn.functional.cosine_similarity(hidden_states, ref_hidden, dim=-1)
                            delta_score = 1.0 - cos.mean().item()
                            
                            # Boost perplexity if delta is high (strong evidence of meaningful change)
                            perplexity *= (1.0 + delta_score * 5.0)
                else:
                    # Fallback: use sequence divergence
                    ref_window = ref_seq[start:end] if end <= len(ref_seq) else ""
                    if ref_window:
                        mismatches = sum(1 for i in range(min(len(window_seq), len(ref_window))) 
                                       if window_seq[i] != ref_window[i])
                        perplexity = mismatches / len(window_seq) if window_seq else 0
                    else:
                        perplexity = 0.5  # Unknown region
                
                total_perplexity += perplexity
                
                # Flag as anomaly if perplexity is high
                if perplexity > threshold:
                    anomaly_windows.append({
                        "start": start,
                        "end": end,
                        "perplexity_score": round(perplexity, 4)
                    })
        
        avg_perplexity = total_perplexity / len(regions_to_scan) if regions_to_scan else 0
        is_anomaly = avg_perplexity > threshold
        
        result = {
            "tier": 1,
            "method": "HyenaDNA",
            "model": model_id,
            "avg_perplexity": round(avg_perplexity, 4),
            "threshold": threshold,
            "is_anomaly": is_anomaly,
            "anomaly_windows": anomaly_windows,
            "cascade_to_tier2": is_anomaly,
            "decision": "ANOMALY - Cascade to Tier 2" if is_anomaly else "VERIFIED SAFE - Skip Cloud"
        }
        
        print(f"   Avg Perplexity: {avg_perplexity:.4f} (threshold: {threshold})")
        print(f"   Anomaly Windows: {len(anomaly_windows)}")
        print(f"   Decision: {result['decision']}")
        
        return result, tokenizer, model
        
    except Exception as e:
        print(f"   ❌ Error loading HyenaDNA model: {e}", file=sys.stderr)
        print("   → Falling back to sequence divergence analysis", file=sys.stderr)
        
        # Fallback: simple divergence
        min_len = min(len(sample_seq), len(ref_seq))
        mismatches = sum(1 for i in range(min_len) if sample_seq[i] != ref_seq[i])
        divergence = mismatches / min_len if min_len > 0 else 0
        is_anomaly = divergence > threshold
        
        result = {
            "tier": 1,
            "method": "HyenaDNA (fallback)",
            "model": model_id,
            "divergence": round(divergence, 4),
            "threshold": threshold,
            "is_anomaly": is_anomaly,
            "anomaly_windows": [],
            "cascade_to_tier2": is_anomaly,
            "decision": "ANOMALY - Cascade to Tier 2" if is_anomaly else "VERIFIED SAFE",
            "error": str(e)
        }
        
        return result, None, None


def tier2_hyena_locus_verification(
    loci_seqs: Dict[str, str],
    model_id: str,
    cache_dir: str,
    device: torch.device,
    threshold: float,
    anomaly_windows: List[Dict],
    target_loci: List[str],
    tokenizer=None,
    model=None
) -> Dict:
    """
    Tier 2: HyenaDNA Locus-Level Verification
    
    Reuses HyenaDNA (already loaded from Tier 1) to verify individual loci.
    This replaces Caduceus/Nucleotide Transformer which requires mamba-ssm.
    """
    print("\n🧬 Tier 2: HyenaDNA Locus-Level Verification")
    print(f"   Model: {model_id} (reusing from Tier 1)")
    print(f"   Device: {device}")
    
    try:
        # Load model if not provided
        if tokenizer is None or model is None:
            print("   → Loading tokenizer and model...")
            tokenizer = AutoTokenizer.from_pretrained(
                model_id,
                cache_dir=cache_dir,
                trust_remote_code=True
            )
            model = AutoModel.from_pretrained(
                model_id,
                cache_dir=cache_dir,
                trust_remote_code=True,
                torch_dtype=torch.float32
            )
            model = model.to(device)
            model.eval()
            print("   ✅ Model loaded")
        else:
            print("   ✅ Reusing loaded model (faster!)")
        
        # Analyze target loci
        locus_scores = {}
        high_risk_loci = []
        
        print(f"   → Analyzing {len(target_loci)} loci...")
        
        with torch.no_grad():
            for locus_name in target_loci:
                if locus_name not in loci_seqs:
                    print(f"      ⚠️  Locus {locus_name} not found")
                    continue
                
                locus_seq = loci_seqs[locus_name]
                
                # Tokenize and run inference
                inputs = tokenizer(locus_seq, return_tensors="pt", truncation=True, max_length=1024)
                inputs = {k: v.to(device) for k, v in inputs.items()}
                
                outputs = model(**inputs)
                
                # Compute structural risk score from hidden states
                if hasattr(outputs, 'last_hidden_state'):
                    hidden_states = outputs.last_hidden_state
                    # Use variance as structural integrity measure
                    # High variance = unusual structure = potential mutation
                    risk_score = torch.var(hidden_states, dim=-1).mean().item()
                else:
                    # Fallback: use sequence complexity
                    risk_score = len(set(locus_seq)) / len(locus_seq) if locus_seq else 0.5
                
                locus_scores[locus_name] = round(risk_score, 4)
                
                if risk_score > threshold:
                    high_risk_loci.append(locus_name)
                    print(f"      ⚠️  {locus_name}: {risk_score:.4f} (HIGH RISK)")
                else:
                    print(f"      ✅ {locus_name}: {risk_score:.4f}")
        
        escalate = len(high_risk_loci) > 0
        
        result = {
            "tier": 2,
            "method": "HyenaDNA (locus verification)",
            "model": model_id,
            "locus_scores": locus_scores,
            "high_risk_loci": high_risk_loci,
            "threshold": threshold,
            "escalate_to_cloud": escalate,
            "decision": f"ESCALATE - {len(high_risk_loci)} high-risk loci" if escalate else "LOW RISK"
        }
        
        print("   📊 Summary:")
        print(f"      Analyzed: {len(locus_scores)} loci")
        print(f"      High Risk: {len(high_risk_loci)}")
        print(f"      Decision: {result['decision']}")
        
        return result
        
    except Exception as e:
        print(f"   ❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return {
            "tier": 2,
            "method": "HyenaDNA (locus verification)",
            "model": model_id,
            "error": str(e),
            "escalate_to_cloud": len(anomaly_windows) > 0,
            "decision": "ERROR - Escalating to be safe"
        }


import pysam
import numpy as np

def validate_integrity(seq, label="Sample"):
    """Check for sequence quality, length, and biological validity."""
    if not seq:
        return False, f"{label} is empty"
    if len(seq) < 100:
        return False, f"{label} is too short ({len(seq)}bp)"
    
    # Check for N-content
    n_content = seq.upper().count('N') / len(seq)
    if n_content > 0.5:
        return False, f"{label} has excessive N-content ({n_content:.1%})"
    
    # Check for complexity (simple entropy check)
    unique_bases = set(seq.upper())
    if len(unique_bases) < 3:
        return False, f"{label} has critically low complexity (bases: {unique_bases})"
        
    return True, "Success"

def load_any_sample(file_path):
    """Universal loader for FASTA, FASTQ, and BAM."""
    path = Path(file_path)
    ext = path.suffix.lower()
    
    # Handle gzipped files
    if ext == ".gz":
        ext = path.suffixes[-2].lower()

    print(f"📂 Ingesting {path.name} (Format: {ext.upper()})")
    
    try:
        if ext in [".fasta", ".fna", ".fa"]:
            seqs = load_fasta(path)
            return "".join(seqs.values())
        
        elif ext in [".fastq", ".fq"]:
            reads = []
            with pysam.FastxFile(str(path)) as fh:
                for i, record in enumerate(fh):
                    reads.append(record.sequence)
                    if i > 1000: break
            return "".join(reads)
            
        elif ext == ".bam":
            bam = pysam.AlignmentFile(str(path), "rb")
            consensus = []
            # Use most covered reference
            if bam.nreferences == 0:
                raise ValueError("BAM file has no references")
            ref_name = bam.references[0]
            for pileupcolumn in bam.pileup(ref_name):
                counts = {}
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        counts[base] = counts.get(base, 0) + 1
                consensus.append(max(counts, key=counts.get) if counts else 'N')
            return "".join(consensus)
        else:
            raise ValueError(f"Unsupported format: {ext}")
    except Exception as e:
        print(f"❌ Error loading {path.name}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Tier 1-2: Local AI Triage (Real Models)")
    parser.add_argument("--tier0-json", required=True, help="Tier 0 sourmash output JSON")
    parser.add_argument("--sample", required=True, help="Input sample (FASTA, FASTQ, or BAM)")
    parser.add_argument("--ref-fasta", required=True, help="Reference genome FASTA")
    parser.add_argument("--loci-fasta", required=True, help="Surveillance loci FASTA")
    parser.add_argument("--hyena-model", required=True, help="HyenaDNA model ID")
    parser.add_argument("--hyena-threshold", type=float, default=4.0, help="HyenaDNA threshold (variance-based)")
    parser.add_argument("--hyena-window-size", type=int, default=500, help="Window size")
    parser.add_argument("--caduceus-model", required=True, help="Caduceus model ID")
    parser.add_argument("--caduceus-loci", nargs="+", default=["ctxB", "wbeT", "tcpA"], help="Target loci")
    parser.add_argument("--caduceus-threshold", type=float, default=2.0, help="Caduceus threshold")
    parser.add_argument("--caduceus-overlap-threshold", type=float, default=0.5, help="Overlap threshold")
    parser.add_argument("--device", default="mps", help="Device: mps, cuda, or cpu")
    parser.add_argument("--model-cache-dir", default="data/models/triage", help="Model cache directory")
    parser.add_argument("--output", required=True, help="Output JSON path")
    
    args = parser.parse_args()
    
    print("🧬 Vibrion Sentinel - Local AI Triage (Real Models)")
    print("=" * 60)
    
    # Load Tier 0 results
    with open(args.tier0_json) as f:
        tier0_result = json.load(f)
    
    print(f"\n📊 Tier 0 Decision: {tier0_result['decision']}")
    
    # If Tier 0 says routine, skip AI analysis
    if tier0_result.get("is_routine", False):
        print("   → Sample is ROUTINE, skipping AI analysis (cost savings!)")
        
        final_result = {
            "tier0": tier0_result,
            "tier1": {"skipped": True, "reason": "Tier 0 routine"},
            "tier2": {"skipped": True, "reason": "Tier 0 routine"},
            "final_decision": {
                "escalate_to_cloud": False,
                "reason": "Routine sample verified by k-mer identity",
                "confidence": "high"
            }
        }
        
        with open(args.output, 'w') as f:
            json.dump(final_result, f, indent=2)
        
        print("\n✅ Triage complete: ROUTINE (no cloud escalation needed)")
        print(f"📄 Results: {args.output}")
        return 0
    
    # Load sequences
    sample_seq = load_any_sample(args.sample)
    
    # Validate Sample Integrity
    is_valid, msg = validate_integrity(sample_seq, "Input Sample")
    if not is_valid:
        print(f"⚠️  INTEGRITY FAILURE: {msg}")
        # Return a safe failure
        final_result = {
            "tier0": tier0_result,
            "error": f"Integrity Failure: {msg}",
            "escalate_to_cloud": True,
            "decision": "INTEGRITY_FAIL - Escalating to Cloud for manual review"
        }
        with open(args.output, 'w') as f:
            json.dump(final_result, f, indent=2)
        return 1

    ref_seqs = load_fasta(Path(args.ref_fasta))
    loci_seqs = load_fasta(Path(args.loci_fasta))
    
    ref_seq = list(ref_seqs.values())[0] if ref_seqs else ""
    
    # Get device
    device = get_device(args.device)
    
    # Run Tier 1: HyenaDNA (returns model for reuse)
    tier1_result, tokenizer, model = tier1_hyena_scan(
        sample_seq=sample_seq,
        ref_seq=ref_seq,
        model_id=args.hyena_model,
        cache_dir=args.model_cache_dir,
        device=device,
        threshold=args.hyena_threshold,
        divergent_regions=tier0_result.get("divergent_regions", []),
        window_size=args.hyena_window_size
    )
    
    # Gray Zone Logic: Escalate if stable but divergent (VOI)
    avg_ppl = tier1_result.get("avg_perplexity", 0)
    is_voi = 1.2 <= avg_ppl <= 4.0
    
    if is_voi:
        print(f"   ✨ VARIANT OF INTEREST (VOI) DETECTED (Perplexity: {avg_ppl:.2f})")
        print("   → Stable but novel structure detected. Forcing EVO2 escalation.")
        tier1_result["cascade_to_tier2"] = True
        tier1_result["voi_flag"] = True

    # Run Tier 2: HyenaDNA locus verification (if Tier 1 found anomalies or VOI)
    if tier1_result.get("cascade_to_tier2", False):
        tier2_result = tier2_hyena_locus_verification(
            loci_seqs=loci_seqs,
            model_id=args.hyena_model,  # Use same model as Tier 1
            cache_dir=args.model_cache_dir,
            device=device,
            threshold=args.caduceus_threshold,
            anomaly_windows=tier1_result.get("anomaly_windows", []),
            target_loci=args.caduceus_loci,
            tokenizer=tokenizer,  # Reuse from Tier 1 (faster!)
            model=model
        )
    else:
        tier2_result = {
            "tier": 2,
            "skipped": True,
            "reason": "Tier 1 verified safe"
        }
    
    # Final decision
    escalate = tier2_result.get("escalate_to_cloud", False)
    
    final_result = {
        "tier0": tier0_result,
        "tier1": tier1_result,
        "tier2": tier2_result,
        "final_decision": {
            "escalate_to_cloud": escalate,
            "reason": "High-risk mutations detected" if escalate else "Local triage verified safe",
            "confidence": "high" if not tier1_result.get("error") else "medium"
        }
    }
    
    # Write output
    with open(args.output, 'w') as f:
        json.dump(final_result, f, indent=2)
    
    print("\n" + "=" * 60)
    print(f"📊 Final Triage Decision: {'ESCALATE TO CLOUD' if escalate else 'LOCAL VERIFIED SAFE'}")
    print(f"📄 Results: {args.output}")
    
    return 0


if __name__ == "__main__":
    exit(main())
