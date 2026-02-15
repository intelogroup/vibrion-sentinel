import torch
import sys
import os

if len(sys.argv) < 2:
    print("Usage: python inspect_ckpt.py <path_to_ckpt>")
    sys.exit(1)

ckpt_path = sys.argv[1]
try:
    print(f"Loading {ckpt_path}...")
    # Add weights_only=False to support legacy/custom objects in checkpoint
    sd = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    if "state_dict" in sd:
        sd = sd["state_dict"]
    
    print("Keys found:")
    for k in sorted(list(sd.keys())): 
        print(k)
        
    print(f"\nTotal keys: {len(sd)}")
    
except Exception as e:
    print(f"Error: {e}")
