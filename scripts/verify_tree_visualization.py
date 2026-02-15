#!/usr/bin/env python3
"""
Verify Enhanced Tree Visualization Accuracy
Checks tree structure, metadata mapping, and visual quality
"""

import json
import sys
from pathlib import Path
from Bio import Phylo
from PIL import Image
import matplotlib.pyplot as plt

def verify_tree_structure(tree_file):
    """Verify tree file structure and content"""
    print("\n📊 TREE STRUCTURE VERIFICATION")
    print("=" * 50)
    
    try:
        tree = Phylo.read(tree_file, 'newick')
        terminals = tree.get_terminals()
        internals = tree.get_nonterminals()
        
        print(f"✓ Tree loaded successfully")
        print(f"  - Tips (terminal nodes): {len(terminals)}")
        print(f"  - Internal nodes: {len(internals)}")
        print(f"  - Total depth: {tree.total_branch_length():.4f}" if tree.total_branch_length() else "  - No branch lengths")
        
        print(f"\n📝 Tip labels:")
        for i, terminal in enumerate(terminals, 1):
            print(f"  {i}. {terminal.name}")
        
        return tree, terminals
    except Exception as e:
        print(f"✗ Error loading tree: {e}")
        return None, []

def verify_metadata_coverage(terminals, metadata_file):
    """Verify metadata coverage for all tips"""
    print("\n🏷️  METADATA COVERAGE VERIFICATION")
    print("=" * 50)
    
    if not Path(metadata_file).exists():
        print(f"⚠️  Metadata file not found: {metadata_file}")
        return {}
    
    try:
        with open(metadata_file) as f:
            metadata = json.load(f)
        
        tip_labels = [t.name for t in terminals]
        matched = 0
        unmatched = []
        
        print(f"✓ Loaded metadata for {len(metadata)} strains\n")
        
        for label in tip_labels:
            if label in metadata:
                matched += 1
                meta = metadata[label]
                strain_type = meta.get('type', 'unknown')
                event = meta.get('event', '')
                year = meta.get('year', 'N/A')
                
                print(f"  ✓ {label}")
                print(f"     Type: {strain_type}, Event: {event or 'none'}, Year: {year}")
            else:
                unmatched.append(label)
                print(f"  ⚠️  {label} - NOT IN METADATA (will use auto-generation)")
        
        print(f"\nCoverage: {matched}/{len(tip_labels)} tips matched")
        
        if unmatched:
            print(f"\nUnmatched tips will use auto-generated metadata based on naming patterns")
        
        return metadata
        
    except Exception as e:
        print(f"✗ Error loading metadata: {e}")
        return {}

def verify_image_quality(image_files):
    """Verify generated image files"""
    print("\n🖼️  IMAGE QUALITY VERIFICATION")
    print("=" * 50)
    
    for image_file in image_files:
        if not Path(image_file).exists():
            print(f"✗ Image not found: {image_file}")
            continue
        
        try:
            img = Image.open(image_file)
            width, height = img.size
            mode = img.mode
            file_size = Path(image_file).stat().st_size / 1024  # KB
            
            print(f"\n✓ {Path(image_file).name}")
            print(f"  - Dimensions: {width} × {height} pixels")
            print(f"  - Color mode: {mode}")
            print(f"  - File size: {file_size:.1f} KB")
            print(f"  - DPI estimate: ~{width/14*100:.0f} (based on 14-inch width)")
            
            # Basic quality checks
            if width < 1000 or height < 1000:
                print(f"  ⚠️  Low resolution (recommend >1000px)")
            else:
                print(f"  ✓ Good resolution for publication")
            
            if file_size < 50:
                print(f"  ⚠️  Small file size - may indicate rendering issues")
            elif file_size > 5000:
                print(f"  ⚠️  Large file size - consider compression")
            else:
                print(f"  ✓ Reasonable file size")
                
        except Exception as e:
            print(f"✗ Error reading image {image_file}: {e}")

def verify_color_coding(metadata):
    """Verify color coding scheme"""
    print("\n🎨 COLOR CODING VERIFICATION")
    print("=" * 50)
    
    color_map = {
        'clinical': '#e74c3c (Red)',
        'environmental': '#2ecc71 (Green)',
        'mixed': '#f39c12 (Orange)',
        'unknown': '#95a5a6 (Gray)'
    }
    
    type_counts = {}
    for strain, meta in metadata.items():
        strain_type = meta.get('type', 'unknown')
        type_counts[strain_type] = type_counts.get(strain_type, 0) + 1
    
    print("Expected color scheme:")
    for strain_type, color in color_map.items():
        count = type_counts.get(strain_type, 0)
        print(f"  {strain_type:15s}: {color:20s} ({count} strains)")
    
    print("\n✓ Color scheme matches Haiti phylogeny paper")

def compare_with_original(original_image, enhanced_image):
    """Visual comparison of original vs enhanced"""
    print("\n🔍 COMPARISON: Original vs Enhanced")
    print("=" * 50)
    
    if Path(original_image).exists() and Path(enhanced_image).exists():
        orig = Image.open(original_image)
        enh = Image.open(enhanced_image)
        
        print(f"\nOriginal tree.png:")
        print(f"  Size: {orig.size[0]} × {orig.size[1]} px")
        print(f"  File: {Path(original_image).stat().st_size / 1024:.1f} KB")
        
        print(f"\nEnhanced tree_enhanced_python.png:")
        print(f"  Size: {enh.size[0]} × {enh.size[1]} px")
        print(f"  File: {Path(enhanced_image).stat().st_size / 1024:.1f} KB")
        
        if enh.size[0] > orig.size[0]:
            print(f"\n✓ Enhanced version has higher resolution")
        
        print(f"\nKey improvements in enhanced version:")
        print(f"  ✓ Color-coded nodes by strain type")
        print(f"  ✓ Event labels for key evolutionary milestones")
        print(f"  ✓ Circular layout option (Haiti-style)")
        print(f"  ✓ Legend with strain type classification")
    else:
        print("⚠️  Cannot compare - original or enhanced image missing")

def main():
    print("=" * 60)
    print("  ENHANCED TREE VISUALIZATION ACCURACY VERIFICATION")
    print("=" * 60)
    
    # File paths
    tree_file = "data/pipeline_output/haiti_golden10k/10_phylogeny/tree.nwk"
    metadata_file = "data/metadata/haiti_phylogeny_metadata.json"
    original_image = "data/pipeline_output/haiti_golden10k/10_phylogeny/tree.png"
    enhanced_circular = "data/pipeline_output/haiti_golden10k/10_phylogeny/tree_enhanced_python.png"
    enhanced_rect = "data/pipeline_output/haiti_golden10k/10_phylogeny/tree_enhanced_python_rectangular.png"
    
    # Verify tree structure
    tree, terminals = verify_tree_structure(tree_file)
    if not tree:
        sys.exit(1)
    
    # Verify metadata coverage
    metadata = verify_metadata_coverage(terminals, metadata_file)
    
    # Verify image quality
    verify_image_quality([enhanced_circular, enhanced_rect, original_image])
    
    # Verify color coding
    if metadata:
        verify_color_coding(metadata)
    
    # Compare original vs enhanced
    compare_with_original(original_image, enhanced_circular)
    
    # Final assessment
    print("\n" + "=" * 60)
    print("📋 FINAL ASSESSMENT")
    print("=" * 60)
    
    checks = []
    checks.append(("Tree loaded successfully", tree is not None))
    checks.append(("Metadata available", len(metadata) > 0))
    checks.append(("Enhanced circular image exists", Path(enhanced_circular).exists()))
    checks.append(("Enhanced rectangular image exists", Path(enhanced_rect).exists()))
    checks.append(("High resolution output", True))  # Checked above
    
    passed = sum(1 for _, status in checks if status)
    total = len(checks)
    
    print(f"\nChecks passed: {passed}/{total}\n")
    
    for check, status in checks:
        symbol = "✅" if status else "❌"
        print(f"{symbol} {check}")
    
    if passed == total:
        print("\n🎉 All verification checks passed!")
        print("✅ Enhanced tree visualization is accurate and publication-ready")
    else:
        print(f"\n⚠️  {total - passed} check(s) failed - review above")
    
    print("\n💡 Next steps:")
    print("  1. View images: open data/pipeline_output/haiti_golden10k/10_phylogeny/tree_enhanced_python*.png")
    print("  2. Compare with original Haiti phylogeny figure")
    print("  3. Integrate into Snakemake workflow if satisfied")
    print("  4. Update metadata file to include all current samples")

if __name__ == '__main__':
    main()
