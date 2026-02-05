#!/usr/bin/env python3
"""
Download and cache HyenaDNA and Caduceus models for local triage.

This script downloads the models to a local cache directory so they don't
need to be re-downloaded every time the pipeline runs.
"""

import argparse
from pathlib import Path
from transformers import AutoModel, AutoTokenizer
import torch

def download_model(model_id: str, cache_dir: Path):
    """Download a model and tokenizer to the cache directory."""
    print(f"\n📥 Downloading {model_id}...")
    print(f"   Cache directory: {cache_dir}")
    
    try:
        # Download tokenizer
        print("   → Downloading tokenizer...")
        tokenizer = AutoTokenizer.from_pretrained(
            model_id,
            cache_dir=str(cache_dir),
            trust_remote_code=True
        )
        print("   ✅ Tokenizer downloaded")
        
        # Download model
        print("   → Downloading model...")
        model = AutoModel.from_pretrained(
            model_id,
            cache_dir=str(cache_dir),
            trust_remote_code=True,
            torch_dtype=torch.float32
        )
        print("   ✅ Model downloaded")
        
        # Get model size
        param_count = sum(p.numel() for p in model.parameters())
        size_mb = param_count * 4 / (1024 * 1024)  # Assuming float32
        print(f"   📊 Model size: {param_count:,} parameters (~{size_mb:.1f} MB)")
        
        return True
        
    except Exception as e:
        print(f"   ❌ Error downloading {model_id}: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Download triage models")
    parser.add_argument(
        "--cache-dir",
        type=str,
        default="data/models/triage",
        help="Directory to cache models"
    )
    parser.add_argument(
        "--hyena-model",
        type=str,
        default="LongSafari/hyenadna-tiny-1k-seqlen-hf",
        help="HyenaDNA model ID (use 'tiny' for testing)"
    )
    parser.add_argument(
        "--caduceus-model",
        type=str,
        default="kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-4",
        help="Caduceus model ID (use smaller variant for testing)"
    )
    
    args = parser.parse_args()
    
    # Create cache directory
    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    print("🧬 Vibrion Sentinel - Triage Model Downloader")
    print("=" * 60)
    print(f"Cache directory: {cache_dir.absolute()}")
    print(f"HyenaDNA model: {args.hyena_model}")
    print(f"Caduceus model: {args.caduceus_model}")
    print("=" * 60)
    
    # Download HyenaDNA
    success_hyena = download_model(args.hyena_model, cache_dir)
    
    # Download Caduceus
    success_caduceus = download_model(args.caduceus_model, cache_dir)
    
    # Summary
    print("\n" + "=" * 60)
    print("📊 Download Summary:")
    print(f"   HyenaDNA: {'✅ Success' if success_hyena else '❌ Failed'}")
    print(f"   Caduceus: {'✅ Success' if success_caduceus else '❌ Failed'}")
    
    if success_hyena and success_caduceus:
        print("\n✅ All models downloaded successfully!")
        print(f"   Models cached in: {cache_dir.absolute()}")
        print("\n💡 Update config.yaml to use these models:")
        print("   triage:")
        print(f"     hyena_model: {args.hyena_model}")
        print(f"     caduceus_model: {args.caduceus_model}")
        print(f"     model_cache_dir: {cache_dir.absolute()}")
    else:
        print("\n⚠️  Some models failed to download. Check errors above.")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
