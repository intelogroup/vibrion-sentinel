"""
Extract unclassified reads from Kraken2 output.
These are reads that didn't match any reference in Kraken2 database.
"""
import sys
import gzip
from Bio import SeqIO

# Snakemake parameters
kraken_out = snakemake.input.kraken_out # noqa: F821
fastq_input = snakemake.input.fastq # noqa: F821
output_fastq = snakemake.output.unclassified # noqa: F821

log_path = snakemake.log[0] if snakemake.log else "/dev/stderr" # noqa: F821
log = open(log_path, "w")

try:
    log.write(f"Parameters: kraken={kraken_out}, fastq={fastq_input}, target={output_fastq}\n")
    # Parse Kraken output and extract unclassified read IDs
    unclassified_ids = set()
    with open(kraken_out, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                classification_status = fields[0]
                read_id = fields[1]
                # 'U' means unclassified
                if classification_status == 'U':
                    unclassified_ids.add(read_id)
    
    log.write(f"Found {len(unclassified_ids)} unclassified reads\n")
    
    # Open input FASTQ (handle both compressed and uncompressed)
    if str(fastq_input).endswith('.gz'):
        input_handle = gzip.open(fastq_input, "rt")
    else:
        input_handle = open(fastq_input, "r")
    
    # Write unclassified reads to output
    if str(output_fastq).endswith('.gz'):
        output_handle = gzip.open(output_fastq, "wt")
    else:
        output_handle = open(output_fastq, "w")
    
    written_count = 0
    for record in SeqIO.parse(input_handle, "fastq"):
        if record.id in unclassified_ids:
            SeqIO.write(record, output_handle, "fastq")
            written_count += 1
    
    output_handle.close()
    input_handle.close()
    
    log.write(f"Wrote {written_count} unclassified reads to {output_fastq}\n")
    log.write("Extraction complete\n")

except Exception as e:
    log.write(f"ERROR: {str(e)}\n")
    sys.exit(1)
finally:
    log.close()
