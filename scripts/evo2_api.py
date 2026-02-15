#!/usr/bin/env python3
"""
EVO2 API interface for strain-level similarity analysis
Uses EVO2-40B generation API for sequence comparison via log-likelihood scoring
"""

import os
import requests
import json
import numpy as np
from typing import List, Dict, Any
from dotenv import load_dotenv

# Load environment variables
load_dotenv(".env.local")

NVIDIA_API_KEY = os.getenv("NVIDIA_API_KEY")
EVO2_GENERATE_URL = "https://health.api.nvidia.com/v1/biology/arc/evo2-40b/generate"

def calculate_sequence_likelihood_evo2(sequence: str, reference: str, max_tokens: int = 100) -> float:
    """
    Calculate how likely a sequence is given a reference using EVO2
    Returns a similarity score between 0 and 1 based on log-likelihood
    
    Higher scores = more similar to reference (same strain/lineage)
    """
    if not NVIDIA_API_KEY:
        raise ValueError("NVIDIA_API_KEY not found in environment")

    headers = {
        "Authorization": f"Bearer {NVIDIA_API_KEY}",
        "Content-Type": "application/json"
    }

    # Use first 500bp for analysis (API token limits)
    analysis_seq = sequence[:500] if len(sequence) > 500 else sequence
    
    payload = {
        "sequence": analysis_seq,
        "num_tokens": max_tokens,
        "temperature": 0.1,
        "top_k": 1,
        "enable_sampled_probs": True
    }

    try:
        response = requests.post(EVO2_GENERATE_URL, json=payload, headers=headers, timeout=60)
        response.raise_for_status()

        data = response.json()
        sampled_probs = data.get("sampled_probs", [])
        
        if not sampled_probs:
            print("⚠️  No probability scores returned")
            return 0.0

        # Average probability = how confident EVO2 is about the sequence
        # Higher confidence = more similar to training data (reference-like sequences)
        avg_prob = sum(sampled_probs) / len(sampled_probs)
        
        return float(avg_prob)

    except requests.exceptions.RequestException as e:
        print(f"❌ EVO2 API request failed: {e}")
        raise
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        raise

def compare_sequences_batch_evo2(query_sequences: List[str], reference_sequence: str) -> List[float]:
    """
    Compare multiple query sequences against a reference
    Returns similarity scores (0-1) for each query sequence
    """
    similarities = []
    
    print(f"Comparing {len(query_sequences)} sequences to reference using EVO2...")
    
    for i, query_seq in enumerate(query_sequences):
        if i % 10 == 0:
            print(f"  Progress: {i}/{len(query_sequences)}")
        
        similarity = calculate_sequence_likelihood_evo2(query_seq, reference_sequence)
        similarities.append(similarity)
    
    print(f"✅ Completed {len(similarities)} comparisons")
    return similarities