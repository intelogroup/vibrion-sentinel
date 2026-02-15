import sys
import os
import subprocess
from pathlib import Path
from Bio import SeqIO

# Mock snakemake for the script
class MockSnakemake:
    def __init__(self, input_ref, input_bed, output_fasta, log_file):
        self.input = type('obj', (object,), {'reference': input_ref, 'bed': input_bed})
        self.output = type('obj', (object,), {'loci_fasta': output_fasta})
        self.log = [log_file]

def run_extraction(ref, bed, out, log):
    import builtins
    builtins.snakemake = MockSnakemake(ref, bed, out, log)
    
    # Import and run the main function from the script
    # We need to add the script's directory to the path
    sys.path.append(os.path.abspath("workflow/scripts"))
    import extract_loci
    extract_loci.main()

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python3 scripts/standalone_extract.py <ref> <bed> <out> <log>")
        sys.exit(1)
    run_extraction(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
