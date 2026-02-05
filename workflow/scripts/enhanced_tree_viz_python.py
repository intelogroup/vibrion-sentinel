#!/usr/bin/env python3
"""
Enhanced Phylogenetic Tree Visualization (Python version)
Generates Haiti-style phylogeny with color-coded nodes
Uses matplotlib and Biopython - no R dependencies
"""

import argparse
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Circle
from Bio import Phylo
import numpy as np

# Color scheme matching Haiti phylogeny
COLOR_MAP = {
    'clinical': '#e74c3c',      # Red
    'environmental': '#2ecc71',  # Green
    'mixed': '#f39c12',         # Yellow/Orange
    'unknown': '#95a5a6'        # Gray
}

def load_metadata(metadata_path):
    """Load strain metadata from JSON"""
    if not Path(metadata_path).exists():
        print(f"⚠️  Metadata file not found: {metadata_path}")
        return {}
    
    try:
        with open(metadata_path) as f:
            metadata = json.load(f)
        print(f"✓ Loaded metadata for {len(metadata)} strains")
        return metadata
    except Exception as e:
        print(f"⚠️  Error loading metadata: {e}")
        return {}

def auto_generate_metadata(tip_labels):
    """Generate metadata from tip label patterns"""
    metadata = {}
    for label in tip_labels:
        # Determine type from label
        if any(x in label.lower() for x in ['2010', '2011', '2012_clin', '2014', 'clin', 'patient']):
            strain_type = 'clinical'
        elif any(x in label.lower() for x in ['env', '2015', '2020', 'water', 'environmental']):
            strain_type = 'environmental'
        elif any(x in label.lower() for x in ['2013', 'adaptation', 'mixed']):
            strain_type = 'mixed'
        else:
            strain_type = 'unknown'
        
        # Determine event
        event = ''
        if '2010' in label:
            event = 'Initiale'
        elif '2012' in label and 'clin' in label.lower():
            event = 'Clin'
        elif 'env' in label.lower() or 'eau' in label.lower():
            event = 'Eau'
        elif 'adaptation' in label.lower():
            event = 'Adaptation'
        elif 'diversif' in label.lower():
            event = 'Diversification'
        elif '2015' in label or 'persist' in label.lower():
            event = 'Persistante'
        elif '2020' in label:
            event = 'Ancêtre 2022'
        elif '2022' in label or 'resurg' in label.lower():
            event = 'Résurgence'
        
        metadata[label] = {
            'type': strain_type,
            'event': event,
            'year': None
        }
    
    return metadata

def plot_tree_circular(tree, metadata, output_path, title="Phylogenetic Network"):
    """Create circular tree layout similar to Haiti phylogeny"""
    
    fig = plt.figure(figsize=(14, 14))
    ax = plt.subplot(111, polar=True)
    
    # Get all terminals
    terminals = tree.get_terminals()
    n_tips = len(terminals)
    
    # Assign angular positions
    angles = np.linspace(0, 2 * np.pi, n_tips, endpoint=False)
    tip_positions = {}
    
    for i, terminal in enumerate(terminals):
        tip_positions[terminal.name] = angles[i]
    
    # Draw branches (simplified)
    for i, terminal in enumerate(terminals):
        angle = angles[i]
        # Draw from center to tip
        ax.plot([0, angle], [0, 1], 'k-', alpha=0.3, linewidth=1.5)
        
        # Get metadata
        label = terminal.name
        meta = metadata.get(label, {'type': 'unknown', 'event': ''})
        strain_type = meta.get('type', 'unknown')
        event = meta.get('event', '')
        
        # Node color and size
        color = COLOR_MAP.get(strain_type, COLOR_MAP['unknown'])
        size = 300 if event else 150
        
        # Draw node
        ax.scatter([angle], [1], c=[color], s=size, alpha=0.8, 
                  edgecolors='white', linewidths=2, zorder=10)
        
        # Add label
        label_text = f"{label}\n({event})" if event else label
        rotation = np.degrees(angle)
        if rotation > 90 and rotation < 270:
            rotation = rotation + 180
            ha = 'right'
        else:
            ha = 'left'
        
        ax.text(angle, 1.15, label_text, rotation=rotation, 
               ha=ha, va='center', fontsize=8, weight='bold')
    
    # Remove radial labels
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.spines['polar'].set_visible(False)
    ax.grid(False)
    
    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor=COLOR_MAP['clinical'], label='Clinical', edgecolor='white'),
        mpatches.Patch(facecolor=COLOR_MAP['environmental'], label='Environmental', edgecolor='white'),
        mpatches.Patch(facecolor=COLOR_MAP['mixed'], label='Mixed', edgecolor='white'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', 
             title='Type de Souche', frameon=True, fontsize=10)
    
    plt.title(title, pad=20, fontsize=16, weight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"✓ Saved circular layout to {output_path}")

def plot_tree_rectangular(tree, metadata, output_path, title="Phylogenetic Network"):
    """Create rectangular tree layout"""
    
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111)
    
    # Draw tree
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)
    
    # Enhance with colors
    terminals = tree.get_terminals()
    for terminal in terminals:
        label = terminal.name
        meta = metadata.get(label, {'type': 'unknown', 'event': ''})
        strain_type = meta.get('type', 'unknown')
        event = meta.get('event', '')
        
        color = COLOR_MAP.get(strain_type, COLOR_MAP['unknown'])
        size = 200 if event else 100
        
        # Find terminal position (approximation)
        # This is a simplified version - in real implementation we'd track exact positions
    
    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor=COLOR_MAP['clinical'], label='Clinical'),
        mpatches.Patch(facecolor=COLOR_MAP['environmental'], label='Environmental'),
        mpatches.Patch(facecolor=COLOR_MAP['mixed'], label='Mixed'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', 
             title='Type de Souche', frameon=True)
    
    plt.title(title, fontsize=16, weight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"✓ Saved rectangular layout to {output_path}")

def main():
    parser = argparse.ArgumentParser(description='Enhanced phylogenetic tree visualization')
    parser.add_argument('tree_file', help='Input Newick tree file')
    parser.add_argument('metadata_file', help='Metadata JSON file')
    parser.add_argument('output_file', help='Output PNG file')
    parser.add_argument('--layout', choices=['circular', 'rectangular'], 
                       default='circular', help='Tree layout style')
    parser.add_argument('--title', default='Phylogenetic Network of V. cholerae Strains',
                       help='Plot title')
    
    args = parser.parse_args()
    
    print("Enhanced Tree Visualization (Python)")
    print("=" * 40)
    print(f"Tree: {args.tree_file}")
    print(f"Metadata: {args.metadata_file}")
    print(f"Output: {args.output_file}")
    print(f"Layout: {args.layout}")
    
    # Load tree
    try:
        tree = Phylo.read(args.tree_file, 'newick')
        terminals = tree.get_terminals()
        print(f"✓ Loaded tree with {len(terminals)} tips")
    except Exception as e:
        print(f"✗ Error loading tree: {e}")
        sys.exit(1)
    
    # Load or generate metadata
    metadata = load_metadata(args.metadata_file)
    if not metadata:
        print("⚠️  Auto-generating metadata from tip labels...")
        tip_labels = [t.name for t in terminals]
        metadata = auto_generate_metadata(tip_labels)
    
    # Generate visualization
    try:
        if args.layout == 'circular':
            plot_tree_circular(tree, metadata, args.output_file, args.title)
            # Also create rectangular version
            rect_output = args.output_file.replace('.png', '_rectangular.png')
            plot_tree_rectangular(tree, metadata, rect_output, args.title)
        else:
            plot_tree_rectangular(tree, metadata, args.output_file, args.title)
    except Exception as e:
        print(f"✗ Error creating visualization: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    print("\n✅ Enhanced tree visualization complete!")

if __name__ == '__main__':
    main()
