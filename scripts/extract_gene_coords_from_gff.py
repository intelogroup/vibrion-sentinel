#!/usr/bin/env python3
"""
Extract gene coordinates from GFF annotation file for the 11 surveillance genes.
"""
import re

# Target genes
GENES = ['wbeT', 'ctxB', 'tcpA', 'gyrA', 'parE', 'scrA', 'scrB', 'scrC', 
         'lip', 'rbmA', 'katB', 'ahpC', 'hapR']

def parse_gff_attributes(attr_string):
    """Parse GFF attributes field"""
    attrs = {}
    for pair in attr_string.split(';'):
        if '=' in pair:
            key, value = pair.split('=', 1)
            attrs[key] = value
    return attrs

def extract_gene_coords(gff_path):
    """Extract coordinates for target genes"""
    gene_coords = {}
    
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chrom, source, feature, start, end, score, strand, phase, attributes = parts
            
            if feature not in ['gene', 'CDS']:
                continue
            
            attrs = parse_gff_attributes(attributes)
            gene_name = attrs.get('gene', attrs.get('Name', ''))
            
            # Check if this is one of our target genes
            for target in GENES:
                if target.lower() == gene_name.lower():
                    if target not in gene_coords or feature == 'gene':
                        gene_coords[target] = {
                            'chrom': chrom,
                            'start': int(start),
                            'end': int(end),
                            'strand': strand,
                            'feature': feature,
                            'locus_tag': attrs.get('locus_tag', 'N/A')
                        }
    
    return gene_coords

if __name__ == "__main__":
    coords = extract_gene_coords("data/references/2010EL-1786.gff")
    
    print("Gene\tChrom\tStart\tEnd\tStrand\tLocus_Tag")
    print("="*80)
    for gene in GENES:
        if gene in coords:
            c = coords[gene]
            print(f"{gene}\t{c['chrom']}\t{c['start']}\t{c['end']}\t{c['strand']}\t{c['locus_tag']}")
        else:
            print(f"{gene}\tNOT FOUND")
    
    print(f"\nFound {len(coords)}/{len(GENES)} genes")
