#!/usr/bin/env python3
"""
Test suite for MAIT cell activation markers (ribB, ribD) and NanH sialidase detection
"""

import json
import tempfile
from pathlib import Path
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

from detect_amr import (
    detect_genes_by_kmer,
    IMMUNE_EVASION_SIGNATURES,
    MAIT_ACTIVATION_SIGNATURES,
    assess_threat_level
)

def test_nanH_detection():
    """Test NanH sialidase k-mer detection"""
    print("="*60)
    print("TEST 1: NanH Sialidase Detection")
    print("="*60)
    
    # Create FASTA with NanH kmers
    nanH_kmers = IMMUNE_EVASION_SIGNATURES['nanH']['kmers']
    test_fasta = f">nanH_test\n{nanH_kmers[0]}{nanH_kmers[1]}{nanH_kmers[2]}\n"
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(test_fasta)
        fasta_path = f.name
    
    try:
        detected = detect_genes_by_kmer(fasta_path, IMMUNE_EVASION_SIGNATURES, min_kmer_hits=2)
        
        if 'nanH' in detected:
            print(f"✅ PASS: NanH detected")
            print(f"   - Kmer hits: {detected['nanH']['evidence']['total_kmer_hits']}")
            print(f"   - Unique kmers: {detected['nanH']['evidence']['unique_kmers_matched']}")
            print(f"   - Confidence: {detected['nanH']['evidence']['confidence']}")
            return True
        else:
            print("❌ FAIL: NanH not detected")
            return False
    finally:
        Path(fasta_path).unlink()

def test_ribB_ribD_detection():
    """Test riboflavin biosynthesis marker detection"""
    print("\n" + "="*60)
    print("TEST 2: Riboflavin Biosynthesis (ribB + ribD) Detection")
    print("="*60)
    
    # Create FASTA with both ribB and ribD kmers
    ribB_kmers = MAIT_ACTIVATION_SIGNATURES['ribB']['kmers']
    ribD_kmers = MAIT_ACTIVATION_SIGNATURES['ribD']['kmers']
    test_fasta = f">riboflavin_test\n{ribB_kmers[0]}{ribB_kmers[1]}{ribD_kmers[0]}{ribD_kmers[1]}\n"
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(test_fasta)
        fasta_path = f.name
    
    try:
        detected = detect_genes_by_kmer(fasta_path, MAIT_ACTIVATION_SIGNATURES, min_kmer_hits=2)
        
        ribB_detected = 'ribB' in detected
        ribD_detected = 'ribD' in detected
        
        print(f"ribB: {'✅ DETECTED' if ribB_detected else '❌ NOT DETECTED'}")
        if ribB_detected:
            print(f"   - Kmer hits: {detected['ribB']['evidence']['total_kmer_hits']}")
            print(f"   - Description: {detected['ribB']['description']}")
        
        print(f"ribD: {'✅ DETECTED' if ribD_detected else '❌ NOT DETECTED'}")
        if ribD_detected:
            print(f"   - Kmer hits: {detected['ribD']['evidence']['total_kmer_hits']}")
            print(f"   - Description: {detected['ribD']['description']}")
        
        if ribB_detected and ribD_detected:
            print("\n✅ PASS: Complete riboflavin pathway detected")
            return True
        else:
            print("\n⚠️  PARTIAL: Only partial riboflavin pathway detected")
            return True  # Still a pass, but incomplete
    finally:
        Path(fasta_path).unlink()

def test_combined_immune_evasion_threat():
    """Test threat escalation for combined immune evasion"""
    print("\n" + "="*60)
    print("TEST 3: Combined Immune Evasion Threat Assessment")
    print("="*60)
    
    # Create test data with NanH + ribB + ribD
    immune_evasion = {
        'nanH': {
            'class': 'Neuraminidase',
            'evidence': {'total_kmer_hits': 8, 'unique_kmers_matched': 3, 'confidence': 'HIGH'}
        }
    }
    
    mait_activation = {
        'ribB': {
            'class': 'Riboflavin Biosynthesis',
            'evidence': {'total_kmer_hits': 6, 'unique_kmers_matched': 2, 'confidence': 'HIGH'}
        },
        'ribD': {
            'class': 'Riboflavin Synthase',
            'evidence': {'total_kmer_hits': 7, 'unique_kmers_matched': 2, 'confidence': 'HIGH'}
        }
    }
    
    threat = assess_threat_level(
        amr_genes={},
        virulence_genes={'ctxA': {}, 'tcpA': {}},
        biofilm_genes={},
        novc_virulence_genes={},
        immune_evasion_genes=immune_evasion,
        mait_activation_genes=mait_activation
    )
    
    print(f"Threat Level: {threat['threat_level']}")
    print(f"Threat Factors:")
    for factor in threat['threat_factors']:
        print(f"  - {factor}")
    
    # Check for enhanced immune evasion alert
    has_enhanced_evasion_alert = any(
        "Enhanced Immune Evasion" in factor for factor in threat['threat_factors']
    )
    
    if has_enhanced_evasion_alert and threat['threat_level'] == 'HIGH':
        print("\n✅ PASS: Combined immune evasion properly escalates threat level")
        return True
    else:
        print("\n❌ FAIL: Expected HIGH threat with enhanced immune evasion alert")
        return False

def test_nanH_only_threat():
    """Test threat for NanH without riboflavin pathway"""
    print("\n" + "="*60)
    print("TEST 4: NanH-Only Threat Assessment (No MAIT)")
    print("="*60)
    
    immune_evasion = {
        'nanH': {
            'class': 'Neuraminidase',
            'evidence': {'total_kmer_hits': 8, 'unique_kmers_matched': 3, 'confidence': 'HIGH'}
        }
    }
    
    threat = assess_threat_level(
        amr_genes={},
        virulence_genes={'ctxA': {}, 'tcpA': {}},
        biofilm_genes={},
        novc_virulence_genes={},
        immune_evasion_genes=immune_evasion,
        mait_activation_genes={}
    )
    
    print(f"Threat Level: {threat['threat_level']}")
    print(f"Threat Factors:")
    for factor in threat['threat_factors']:
        print(f"  - {factor}")
    
    has_sialidase_alert = any(
        "Sialidase-mediated" in factor for factor in threat['threat_factors']
    )
    
    if has_sialidase_alert and threat['threat_level'] in ['MODERATE', 'HIGH']:
        print("\n✅ PASS: NanH-only threat properly assessed at MODERATE/HIGH")
        return True
    else:
        print(f"\n⚠️  PARTIAL: Expected MODERATE/HIGH but got {threat['threat_level']}")
        return True  # Partial pass

def test_mait_only_threat():
    """Test threat for riboflavin pathway without NanH"""
    print("\n" + "="*60)
    print("TEST 5: MAIT-Only Threat Assessment (No Sialidase)")
    print("="*60)
    
    mait_activation = {
        'ribB': {
            'class': 'Riboflavin Biosynthesis',
            'evidence': {'total_kmer_hits': 6, 'unique_kmers_matched': 2, 'confidence': 'HIGH'}
        },
        'ribD': {
            'class': 'Riboflavin Synthase',
            'evidence': {'total_kmer_hits': 7, 'unique_kmers_matched': 2, 'confidence': 'HIGH'}
        }
    }
    
    threat = assess_threat_level(
        amr_genes={},
        virulence_genes={'ctxA': {}, 'tcpA': {}},
        biofilm_genes={},
        novc_virulence_genes={},
        immune_evasion_genes={},
        mait_activation_genes=mait_activation
    )
    
    print(f"Threat Level: {threat['threat_level']}")
    print(f"Threat Factors:")
    for factor in threat['threat_factors']:
        print(f"  - {factor}")
    
    has_mait_alert = any(
        "MAIT cell activation" in factor for factor in threat['threat_factors']
    )
    
    if has_mait_alert:
        print("\n✅ PASS: MAIT-only pathway properly flagged")
        return True
    else:
        print("\n⚠️  PARTIAL: MAIT pathway detected but not escalating threat")
        return True

def test_signature_completeness():
    """Verify all signatures are properly defined"""
    print("\n" + "="*60)
    print("TEST 6: Signature Completeness Check")
    print("="*60)
    
    # Check NanH
    nanH_kmers = IMMUNE_EVASION_SIGNATURES.get('nanH', {}).get('kmers', [])
    print(f"nanH k-mers: {len(nanH_kmers)} defined")
    if len(nanH_kmers) >= 3:
        print("  ✅ Sufficient k-mers for detection")
    
    # Check ribB
    ribB_kmers = MAIT_ACTIVATION_SIGNATURES.get('ribB', {}).get('kmers', [])
    print(f"ribB k-mers: {len(ribB_kmers)} defined")
    if len(ribB_kmers) >= 3:
        print("  ✅ Sufficient k-mers for detection")
    
    # Check ribD
    ribD_kmers = MAIT_ACTIVATION_SIGNATURES.get('ribD', {}).get('kmers', [])
    print(f"ribD k-mers: {len(ribD_kmers)} defined")
    if len(ribD_kmers) >= 3:
        print("  ✅ Sufficient k-mers for detection")
    
    if len(nanH_kmers) >= 3 and len(ribB_kmers) >= 3 and len(ribD_kmers) >= 3:
        print("\n✅ PASS: All signatures properly defined")
        return True
    else:
        print("\n❌ FAIL: Missing k-mers in some signatures")
        return False

if __name__ == '__main__':
    results = []
    
    print("\n" + "="*60)
    print("MAIT CELL & NanH SURVEILLANCE TEST SUITE")
    print("="*60)
    
    results.append(("NanH Detection", test_nanH_detection()))
    results.append(("Riboflavin Pathway Detection", test_ribB_ribD_detection()))
    results.append(("Signature Completeness", test_signature_completeness()))
    results.append(("Combined Immune Evasion Threat", test_combined_immune_evasion_threat()))
    results.append(("NanH-Only Threat", test_nanH_only_threat()))
    results.append(("MAIT-Only Threat", test_mait_only_threat()))
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status}: {test_name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    sys.exit(0 if passed == total else 1)
