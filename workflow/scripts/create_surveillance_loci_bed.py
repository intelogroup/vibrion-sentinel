#!/usr/bin/env python3
"""
Create BED file with absolute genomic coordinates for surveillance loci.
Converts relative positions (0.0-1.0) to actual chromosome coordinates.
"""
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from backend.core.loci import (
    SURVEILLANCE_LOCI_FULL,
    CHR1_LENGTH
)

# Window size for each locus (500bp on each side of center)
WINDOW_SIZE = 1000


def create_bed_file():
    """Generate BED file with absolute coordinates"""
    
    output_path = Path(__file__).parent.parent.parent / "data" / "references" / "surveillance_loci.bed"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print("📝 Creating surveillance loci BED file")
    print(f"   Output: {output_path}")
    
    with open(output_path, "w") as bed:
        # BED format: chromosome, start, end, name
        bed.write("# Surveillance loci for Vibrio cholerae 2010EL-1786 reference\n")
        bed.write("# Format: chromosome start end name\n")
        
        for locus_name, relative_pos in SURVEILLANCE_LOCI_FULL:
            # All loci are on chromosome 1 for now
            # Calculate absolute position
            center_pos = int(relative_pos * CHR1_LENGTH)
            
            # Create 1000bp window centered on position
            start = max(0, center_pos - WINDOW_SIZE // 2)
            end = min(CHR1_LENGTH, center_pos + WINDOW_SIZE // 2)
            
            # Write BED entry
            bed.write(f"CP003069.1\t{start}\t{end}\t{locus_name}\n")
            
            print(f"   ✓ {locus_name}: {start:,} - {end:,} bp (center: {center_pos:,})")
    
    print(f"\n✅ Created BED file with {len(SURVEILLANCE_LOCI_FULL)} loci")
    return output_path


if __name__ == "__main__":
    create_bed_file()
