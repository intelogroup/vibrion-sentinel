"""
Merge classified Vibrio reads with NT-500M rescued reads.
Ensures complete coverage of all Vibrio DNA in sample.
"""
import sys
import gzip
from Bio import SeqIO

classified_fastq = snakemake.input.classified # noqa: F821
rescued_fastq = snakemake.input.rescued # noqa: F821
output_fastq = snakemake.output.merged # noqa: F821

log = open(snakemake.log[0], "w") # noqa: F821

try:
    # Open input files
    if str(classified_fastq).endswith('.gz'):
        classified_handle = gzip.open(classified_fastq, "rt")
    else:
        classified_handle = open(classified_fastq, "r")
    
    if str(rescued_fastq).endswith('.gz'):
        rescued_handle = gzip.open(rescued_fastq, "rt")
    else:
        rescued_handle = open(rescued_fastq, "r")
    
    # Open output file
    if str(output_fastq).endswith('.gz'):
        output_handle = gzip.open(output_fastq, "wt")
    else:
        output_handle = open(output_fastq, "w")
    
    # Write classified reads
    classified_count = 0
    for record in SeqIO.parse(classified_handle, "fastq"):
        SeqIO.write(record, output_handle, "fastq")
        classified_count += 1
    
    classified_handle.close()
    
    # Write rescued reads
    rescued_count = 0
    for record in SeqIO.parse(rescued_handle, "fastq"):
        SeqIO.write(record, output_handle, "fastq")
        rescued_count += 1
    
    rescued_handle.close()
    output_handle.close()
    
    total_count = classified_count + rescued_count
    log.write(f"Merged {classified_count} classified + {rescued_count} rescued = {total_count} total reads\n")
    log.write(f"Output: {output_fastq}\n")

except Exception as e:
    log.write(f"ERROR: {str(e)}\n")
    sys.exit(1)
finally:
    log.close()
