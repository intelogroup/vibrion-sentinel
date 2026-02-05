#!/usr/bin/env python3
"""
Vibrion Sentinel: Prep Phylogeny
Combines sample consensus with reference library for MSA.
"""

import argparse
from pathlib import Path
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--consensus", required=True, help="Sample consensus FASTA")
    parser.add_argument("--sample-id", required=True, help="Sample ID")
    parser.add_argument("--output", required=True, help="Output Combined FASTA")
    args = parser.parse_args()
    
    # 1. Load Sample
    # Note: Consensus usually has one contig or is the Ref with SNPs applied
    sample_recs = list(SeqIO.parse(args.consensus, "fasta"))
    # Rename specifically for the tree
    for i, rec in enumerate(sample_recs):
        rec.id = args.sample_id
        rec.description = f"Sample_{args.sample_id}"
        
    # 2. Load References
    # We use the 'production_references' produced/downloaded earlier
    ref_dir = Path("data/production_references")
    
    references = [
        ("Haiti_2010", ref_dir / "Haiti_2010_Ref.fasta"), # Or 2010EL-1786.fasta
        ("Bengal_1993", ref_dir / "Bengal_1993_Combined.fasta"),
        ("Classical_569B", ref_dir / "Classical_569B_Combined.fasta"),
        ("Environmental_NOVC", ref_dir / "Environmental_NOVC_Contig1.fasta")
    ]
    
    final_records = []
    final_records.extend(sample_recs)
    
    for name, path in references:
        if path.exists():
            # For multi-contig references, we might just take the largest or concatenate
            # For a quick tree, taking the largest component is often safest for alignment if not doing whole-genome alignment
            # But let's try to just dump them. Mafft usually handles it or we treat as unaligned.
            # Actually, `augur align` or `mafft` prefers 1 sequence per taxon.
            # We should probably concatenate if multi-fasta.
            
            recs = list(SeqIO.parse(path, "fasta"))
            full_seq = "".join([str(r.seq) for r in recs])
            
            # Simple record creation
            from Bio.SeqRecord import SeqRecord
            from Bio.Seq import Seq
            
            new_rec = SeqRecord(
                Seq(full_seq),
                id=name,
                description=f"Reference_Archetype_{name}"
            )
            final_records.append(new_rec)
        else:
            # Try alternative path (standard references)
            alt_path = Path("data/references") / "2010EL-1786.fasta"
            if name == "Haiti_2010" and alt_path.exists():
                 recs = list(SeqIO.parse(alt_path, "fasta"))
                 full_seq = "".join([str(r.seq) for r in recs])
                 new_rec = SeqRecord(Seq(full_seq), id=name, description=f"Reference_Archetype_{name}")
                 final_records.append(new_rec)
            else:
                print(f"Warning: Reference {name} not found at {path}")

    # 3. Write Output
    SeqIO.write(final_records, args.output, "fasta")
    print(f"Prepared {len(final_records)} sequences for phylogeny.")

if __name__ == "__main__":
    main()
