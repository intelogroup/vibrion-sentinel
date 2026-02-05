#!/usr/bin/env python3
"""
Haiti-Specific Phylogenetic Tree Builder (2010-2022)
Creates a focused phylogeny of Haiti cholera outbreak evolution
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess
import json

def extract_core_loci(fasta_file, bed_file, strain_name):
    """Extract core surveillance loci from a genome"""
    try:
        # Read genome
        genome = SeqIO.read(fasta_file, "fasta")
        
        # Read BED coordinates
        loci_sequences = []
        with open(bed_file) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                    
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                locus_name = parts[3] if len(parts) > 3 else f"{chrom}_{start}_{end}"
                
                # Extract sequence
                locus_seq = str(genome.seq[start:end])
                loci_sequences.append(locus_seq)
        
        # Concatenate all loci
        concatenated = "".join(loci_sequences)
        return SeqRecord(Seq(concatenated), id=strain_name, description=f"Core loci from {fasta_file}")
        
    except Exception as e:
        print(f"Warning: Could not extract loci from {fasta_file}: {e}")
        return None

def create_haiti_alignment(output_dir):
    """Create alignment of Haiti strains 2010-2022"""
    
    # Use pre-extracted loci files
    haiti_loci = {
        "Haiti_2010_Ancestor": "data/references/2010EL-1786.fasta",  # Will extract
        "Haiti_2022_Resurgent": "data/references/Haiti_2022_Resurgence_loci_v2.fasta",
        "Bengal_1993": "data/references/bengal_1993_loci.fasta",
        "Classical_569B": "data/references/classical_569b_loci.fasta"
    }
    
    bed_file = "data/references/surveillance_loci.bed"
    alignment_records = []
    
    print("Loading Haiti strain loci...")
    for strain_name, fasta_path in haiti_loci.items():
        if not Path(fasta_path).exists():
            print(f"⚠️  Warning: {fasta_path} not found, skipping {strain_name}")
            continue
        
        # If it's a full genome, extract loci; otherwise load directly
        if "loci" in fasta_path:
            try:
                # Concatenate multi-locus file
                loci_records = list(SeqIO.parse(fasta_path, "fasta"))
                if loci_records:
                    combined_seq = "".join(str(rec.seq) for rec in loci_records)
                    record = SeqRecord(Seq(combined_seq), id=strain_name, 
                                     description=f"Core loci from {fasta_path}")
                    alignment_records.append(record)
                    print(f"  ✓ {strain_name}: {len(record.seq)} bp ({len(loci_records)} loci)")
            except Exception as e:
                print(f"  ✗ Error loading {strain_name}: {e}")
        else:
            # Extract loci from full genome
            record = extract_core_loci(fasta_path, bed_file, strain_name)
            if record:
                alignment_records.append(record)
                print(f"  ✓ {strain_name}: {len(record.seq)} bp (extracted)")
    
    if len(alignment_records) < 2:
        print("Error: Need at least 2 strains for phylogeny")
        return None
    
    # Add test sample if available
    sample_file = output_dir / "sample_core.fasta"
    if sample_file.exists():
        try:
            sample_rec = SeqIO.read(sample_file, "fasta")
            alignment_records.append(sample_rec)
            print(f"  ✓ {sample_rec.id}: {len(sample_rec.seq)} bp (test sample)")
        except:
            pass
    
    # Create alignment file
    alignment_file = output_dir / "haiti_core_alignment.fasta"
    SeqIO.write(alignment_records, alignment_file, "fasta")
    print(f"\n✓ Created unaligned sequences: {alignment_file}")
    print(f"  Strains: {len(alignment_records)}")
    
    # Align with MAFFT
    print("\nAligning sequences with MAFFT...")
    aligned_file = output_dir / "haiti_core_aligned.fasta"
    try:
        subprocess.run([
            "mafft", "--auto", "--thread", "4",
            str(alignment_file)
        ], check=True, stdout=open(aligned_file, 'w'), stderr=subprocess.DEVNULL)
        print(f"  ✓ Alignment complete: {aligned_file}")
        return aligned_file, alignment_records
    except Exception as e:
        print(f"  ⚠️  MAFFT failed: {e}")
        print("  Using unaligned sequences (may cause issues)")
        return alignment_file, alignment_records

def build_tree_iqtree(alignment_file, output_prefix):
    """Build ML tree using IQ-TREE"""
    
    # Check for iqtree
    iqtree_path = None
    for cmd in ['iqtree2', 'iqtree']:
        try:
            subprocess.run([cmd, '-v'], capture_output=True, check=True)
            iqtree_path = cmd
            break
        except:
            continue
    
    if not iqtree_path:
        print("⚠️  IQ-TREE not found, using distance-based method")
        return None
    
    print(f"\nBuilding ML tree with {iqtree_path}...")
    cmd = [
        iqtree_path,
        '-s', str(alignment_file),
        '-m', 'GTR+G',
        '-nt', '4',
        '-pre', str(output_prefix),
        '-redo'
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        tree_file = Path(f"{output_prefix}.treefile")
        if tree_file.exists():
            print("✓ ML tree built successfully")
            return tree_file
    except Exception as e:
        print(f"⚠️  IQ-TREE failed: {e}")
        return None

def build_tree_upgma(alignment_file, output_file):
    """Build UPGMA tree as fallback"""
    
    print("\nBuilding UPGMA tree...")
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        
        # Calculate distance matrix
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        
        # Build tree
        constructor = DistanceTreeConstructor(calculator)
        tree = constructor.upgma(dm)
        
        # Write tree
        Phylo.write(tree, output_file, "newick")
        print(f"✓ UPGMA tree saved: {output_file}")
        return output_file
        
    except Exception as e:
        print(f"✗ Error building UPGMA tree: {e}")
        return None

def create_haiti_metadata():
    """Create Haiti-specific metadata"""
    
    metadata = {
        "Haiti_2010_Ancestor": {
            "type": "clinical",
            "event": "Initiale (2010)",
            "year": 2010,
            "location": "Artibonite Valley",
            "source": "Outbreak index case region",
            "description": "2010EL-1786 - Original outbreak ancestor"
        },
        "Haiti_2022_Resurgent": {
            "type": "clinical",
            "event": "Résurgence (2022)",
            "year": 2022,
            "location": "Port-au-Prince",
            "source": "Post-silence resurgence",
            "description": "2022 outbreak resurgence after 3-year silence"
        },
        "haiti_golden10k": {
            "type": "clinical",
            "event": "",
            "year": 2024,
            "location": "Unknown",
            "source": "Test sample",
            "description": "Current pipeline test sample"
        }
    }
    
    return metadata

def main():
    parser = argparse.ArgumentParser(description='Build Haiti-specific phylogenetic tree (2010-2022)')
    parser.add_argument('--sample-dir', required=True, help='Sample output directory')
    parser.add_argument('--output', required=True, help='Output tree file')
    parser.add_argument('--method', choices=['iqtree', 'upgma', 'auto'], default='auto',
                       help='Tree building method')
    
    args = parser.parse_args()
    
    output_dir = Path(args.sample_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("  HAITI CHOLERA PHYLOGENY BUILDER (2010-2022)")
    print("=" * 60)
    
    # Create alignment
    alignment_file, records = create_haiti_alignment(output_dir)
    if not alignment_file:
        sys.exit(1)
    
    # Build tree
    tree_file = None
    
    if args.method in ['iqtree', 'auto']:
        tree_prefix = output_dir / "haiti_tree"
        tree_file = build_tree_iqtree(alignment_file, tree_prefix)
    
    if not tree_file and args.method in ['upgma', 'auto']:
        tree_file = build_tree_upgma(alignment_file, output_dir / "haiti_tree_upgma.nwk")
    
    if not tree_file:
        print("✗ Failed to build tree")
        sys.exit(1)
    
    # Copy to final output
    import shutil
    shutil.copy(tree_file, args.output)
    print(f"\n✓ Final tree saved: {args.output}")
    
    # Save metadata
    metadata = create_haiti_metadata()
    metadata_file = output_dir / "haiti_metadata.json"
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"✓ Metadata saved: {metadata_file}")
    
    # Print tree
    print("\n" + "=" * 60)
    print("TREE STRUCTURE")
    print("=" * 60)
    tree = Phylo.read(args.output, "newick")
    Phylo.draw_ascii(tree)
    
    print(f"\n✅ Haiti phylogeny complete!")
    print(f"   Strains: {len(records)}")
    print(f"   Tree: {args.output}")
    print(f"   Metadata: {metadata_file}")

if __name__ == '__main__':
    main()
