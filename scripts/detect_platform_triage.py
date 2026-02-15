#!/usr/bin/env python3
"""
Vibrion Sentinel: Platform Triage
Detects compute resources and configures pipeline physics.
Implements "Graceful Degradation" for field operations.
"""

import shutil
import json
import sys
import subprocess

def check_nvidia_smi():
    """Check for NVIDIA GPU via nvidia-smi."""
    if shutil.which("nvidia-smi") is None:
        return False, "Not Found"
    
    try:
        # Check if we can actually query it (detects driver issues)
        subprocess.check_output(["nvidia-smi", "-L"], stderr=subprocess.STDOUT)
        return True, "Active"
    except subprocess.CalledProcessError:
        return False, "Driver Error"

def main():
    # default fallback
    config = {
        "platform": "CPU",
        "medaka_device": "cpu",
        "mmseqs_threads": 4, # Conservative default
        "guppy_device": "cpu", # If we were basecalling
        "warning": None
    }
    
    # 1. Check GPU
    has_gpu, status = check_nvidia_smi()
    
    if has_gpu:
        config["platform"] = "CUDA"
        config["medaka_device"] = "cuda:0"
        config["warning"] = None
        # Optimization: If GPU exists, we can probably spare more CPU threads for IO
        config["mmseqs_threads"] = 8 
    else:
        config["platform"] = "CPU_FALLBACK"
        config["medaka_device"] = "cpu"
        config["warning"] = f"GPU Unavailable ({status}). Running in Bunker Mode (Slow)."

    # 2. Check RAM (Optional - adjust threads if low RAM?)
    # ...

    # Output JSON for Snakemake or logs
    print(json.dumps(config, indent=2))
    
    if config["warning"]:
        print(f"⚠️  {config['warning']}", file=sys.stderr)

if __name__ == "__main__":
    main()
