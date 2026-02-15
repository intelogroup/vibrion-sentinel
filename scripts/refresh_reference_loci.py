import sys
import subprocess
from pathlib import Path

def extract_loci(reference_path, bed_path, output_fasta):
    reference_path = Path(reference_path)
    bed_path = Path(bed_path)
    
    # Index reference
    subprocess.run(['samtools', 'faidx', str(reference_path)], check=True)
    
    with open(output_fasta, 'w') as out:
        with open(bed_path) as bed:
            for line in bed:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 4: continue
                
                chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
                
                # samtools faidx uses 1-based inclusive
                # BED is 0-based exclusive-end
                # So start+1 to end is correct for 1-based
                coord = f"{chrom}:{start+1}-{end}"
                
                try:
                    res = subprocess.run(['samtools', 'faidx', str(reference_path), coord], capture_output=True, text=True, check=True)
                    seq = "".join(res.stdout.split('\n')[1:])
                    out.write(f">{name}\n{seq}\n")
                    print(f"Extracted {name}")
                except Exception as e:
                    print(f"Failed to extract {name}: {e}")

if __name__ == "__main__":
    extract_loci("data/references/2010EL-1786.fasta", "data/references/surveillance_loci.bed", "data/references/reference_loci.fasta")

