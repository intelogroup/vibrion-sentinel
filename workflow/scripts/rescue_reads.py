"""
Calculate cosine similarity between unclassified read embeddings
and Vibrio reference embedding. Rescue reads >= threshold.
"""
import sys
import json
import gzip
import numpy as np
from Bio import SeqIO
from sklearn.metrics.pairwise import cosine_similarity

unclassified_fastq = snakemake.input.unclassified # noqa: F821
embeddings_path = snakemake.input.embeddings # noqa: F821
read_ids_path = snakemake.input.read_ids # noqa: F821
ref_embedding_path = snakemake.input.ref_embedding # noqa: F821
output_fastq = snakemake.output.rescued # noqa: F821
output_stats = snakemake.output.stats # noqa: F821
threshold = snakemake.params.similarity_threshold # noqa: F821

log = open(snakemake.log[0], "w") # noqa: F821

try:
    # Load embeddings
    read_embeddings = np.load(embeddings_path)
    ref_embedding = np.load(ref_embedding_path).reshape(1, -1)
    
    # Load read IDs from JSON file passed as input
    with open(read_ids_path, 'r') as f:
        read_ids = json.load(f)
    
    if len(read_embeddings) == 0:
        log.write("No unclassified reads to rescue\n")
        # Create empty output
        if str(output_fastq).endswith('.gz'):
            gzip.open(output_fastq, "wt").close()
        else:
            open(output_fastq, "w").close()
        
        stats = {
            "total_reads": 0,
            "rescued_reads": 0,
            "similarity_threshold": threshold,
            "min_similarity": None,
            "max_similarity": None
        }
        with open(output_stats, 'w') as f:
            json.dump(stats, f, indent=2)
        log.close()
        sys.exit(0)
    
    # Calculate cosine similarities
    similarities = cosine_similarity(read_embeddings, ref_embedding).flatten()
    
    # Find reads to rescue
    rescue_mask = similarities >= threshold
    rescued_ids = set([read_ids[i] for i in range(len(read_ids)) if rescue_mask[i]])
    
    log.write(f"Rescued {len(rescued_ids)} / {len(read_ids)} unclassified reads (threshold={threshold})\n")
    log.write(f"Similarity range: {similarities.min():.4f} - {similarities.max():.4f}\n")
    log.write(f"Mean similarity: {similarities.mean():.4f}\n")
    log.write(f"Median similarity: {np.median(similarities):.4f}\n")
    
    # Write rescued reads to FASTQ
    if str(unclassified_fastq).endswith('.gz'):
        input_handle = gzip.open(unclassified_fastq, "rt")
    else:
        input_handle = open(unclassified_fastq, "r")
    
    if str(output_fastq).endswith('.gz'):
        output_handle = gzip.open(output_fastq, "wt")
    else:
        output_handle = open(output_fastq, "w")
    
    written_count = 0
    for record in SeqIO.parse(input_handle, "fastq"):
        if record.id in rescued_ids:
            SeqIO.write(record, output_handle, "fastq")
            written_count += 1
    
    output_handle.close()
    input_handle.close()
    
    log.write(f"Wrote {written_count} rescued reads to {output_fastq}\n")
    
    # Write statistics
    stats = {
        "total_unclassified_reads": len(read_ids),
        "rescued_reads": len(rescued_ids),
        "rescue_percentage": round((len(rescued_ids) / len(read_ids) * 100), 2),
        "similarity_threshold": threshold,
        "min_similarity": float(similarities.min()),
        "max_similarity": float(similarities.max()),
        "mean_similarity": float(similarities.mean()),
        "median_similarity": float(np.median(similarities)),
        "std_similarity": float(similarities.std())
    }
    
    with open(output_stats, 'w') as f:
        json.dump(stats, f, indent=2)
    
    log.write("Rescue complete\n")

except Exception as e:
    log.write(f"ERROR: {str(e)}\n")
    import traceback
    log.write(traceback.format_exc())
    sys.exit(1)
finally:
    log.close()
