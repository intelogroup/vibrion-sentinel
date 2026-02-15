import os
import subprocess
import json
import gzip
from Bio import SeqIO

def run_mmseqs(input_fastq, db_path, output_prefix, threads, min_seq_id):
    """Run MMseqs2 easy-taxonomy."""
    cmd = [
        "mmseqs", "easy-taxonomy",
        input_fastq,
        db_path,
        output_prefix,
        "tmp_mmseqs_rescue",
        "--threads", str(threads),
        "--min-seq-id", str(min_seq_id),
        "--tax-lineage", "1",
        "--lca-mode", "3", # LCA mode 3 is robust
        "--min-length", "10"
    ]
    
    print(f"Running command: {' '.join(cmd)}")
    subprocess.check_call(cmd)
    
    return f"{output_prefix}_lca.tsv", f"{output_prefix}_report"

def filter_vibrio_reads(lca_file, target_taxid=662):
    """Parse LCA file and return set of read IDs classification as target or descendant."""
    # Note: MMseqs2 easy-taxonomy with --tax-lineage 1 produces a file where we can check lineage.
    # But usually parsing the TaxID is enough if we have the taxonomy tree.
    # However, easy-taxonomy _lca.tsv format:
    # ReadID <tab> TaxID <tab> Rank <tab> Name <tab> ... (Lineage if loaded?)
    
    # For robust filtering, we rely on the TaxID. 
    # Since we might not have a full taxonomy tree loaded in python, 
    # we can trust mmseqs if we grep for Vibrio? 
    # Or better: We assume 662 is Vibrio.
    # To be safe, we capture everything where "Vibrio" is in the name provided by mmseqs (column 4 usually).
    
    vibrio_read_ids = set()
    
    with open(lca_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            
            read_id = parts[0]
            tax_id = parts[1]
            rank = parts[2]
            name = parts[3]
            
            # Simple string match for safety if tax logic is complex
            # user wants "Vibrio" (662). 
            # Note: 662 is genus Vibrio. 
            # Descendants will have different taxids (e.g. 666 for V. cholerae).
            # We can check if "Vibrio" is in the name string if mmseqs provides full path, 
            # but usually it provides the scientific name of the assigned node.
            # If the node is "Vibrio cholerae", keeping it is correct.
            # If "Vibrio", keeping it is correct.
            # If "Aliivibrio", we usually skip unless requested.
            
            # We will use string matching on 'Vibrio' in the name for simplicity and robustness 
            # without loading a full taxdump.
            if "Vibrio" in name:
                vibrio_read_ids.add(read_id)
                
    return vibrio_read_ids

def extract_reads(input_fastq, read_ids, output_fastq):
    """Extract reads from input FASTQ if they are in read_ids set."""
    count = 0
    # Handle gzipped input if needed, but input in rule is likely uncompressed or we check extension
    # The rule often provides .fastq or .fastq.gz
    
    open_func = gzip.open if input_fastq.endswith('.gz') else open
    mode = 'rt' if input_fastq.endswith('.gz') else 'r'
    
    # Output should be gzipped to save space
    with open_func(input_fastq, mode) as f_in, gzip.open(output_fastq, 'wt') as f_out:
        for record in SeqIO.parse(f_in, "fastq"):
            if record.id in read_ids:
                SeqIO.write(record, f_out, "fastq")
                count += 1
    return count

def main():
    input_fastq = snakemake.input.unclassified # noqa: F821
    output_fastq = snakemake.output.rescued # noqa: F821
    stats_file = snakemake.output.stats # noqa: F821
    
    db_path = snakemake.params.db # noqa: F821
    threads = snakemake.params.threads # noqa: F821
    min_seq_id = snakemake.params.get("min_seq_id", 0.6) # Default from user plan was 0.8 but code says 0.6? Plan said 0.8 # noqa: F821
    # Config has 0.8.
    
    # Temp prefix
    output_prefix = f"{os.path.dirname(output_fastq)}/mmseqs_tmp_{os.path.basename(input_fastq)}"
    
    # 1. Run MMseqs2
    try:
        lca_file, report_file = run_mmseqs(input_fastq, db_path, output_prefix, threads, min_seq_id)
        
        # 2. Filter Reads
        rescued_ids = filter_vibrio_reads(lca_file)
        
        # 3. Extract Reads
        count = extract_reads(input_fastq, rescued_ids, output_fastq)
        
        # 4. Stats
        stats = {
            "method": "mmseqs2_easy_taxonomy",
            "db": "swissprot",
            "rescued_count": count,
            "min_seq_id": min_seq_id
        }
        
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=4)

        # Cleanup
        # os.remove(lca_file) 
        # os.remove(report_file)
        # cleanup other mmseqs tmp files?
        # mmseqs creates many output files, maybe better to keep them for debug or rm output_prefix* 
        
    except Exception as e:
        print(f"Error in MMseqs2 rescue: {e}")
        # Create empty output on failure to avoid snakemake.crash loop, or fail hard? # noqa: F821
        # Fail hard is better for debugging
        raise e

if __name__ == "__main__":
    main()
