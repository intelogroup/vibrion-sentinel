
import os
import sys
import json
import gzip
import argparse
import asyncio
import aiohttp
from pathlib import Path
from Bio import SeqIO
from dotenv import load_dotenv

# Load environment variables from project root
root_path = Path(__file__).parent.parent.parent
load_dotenv(root_path / ".env.local")

NVIDIA_API_KEY = os.getenv("NVIDIA_API_KEY")
EVO2_GENERATE_URL = "https://health.api.nvidia.com/v1/biology/arc/evo2-40b/generate"


def _open_text_maybe_gzip(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _get_io_from_snakemake_or_cli():
    """Supports both Snakemake script execution and standalone CLI usage."""
    if "snakemake" in globals():
        return {
            "unclassified_fastq": snakemake.input.unclassified,
            "output_fastq": snakemake.output.rescued,
            "output_stats": snakemake.output.stats,
            "threshold": getattr(snakemake.params, "confidence_threshold", 0.6),
            "max_reads_to_process": getattr(snakemake.params, "max_reads", 5000),
            "evo2_url": getattr(snakemake.params, "evo2_url", "https://health.api.nvidia.com/v1/biology/arc/evo2-7b/generate"),
            "log_path": snakemake.log[0] if snakemake.log else None,
        }

    parser = argparse.ArgumentParser(description="Evo2 API read rescue")
    parser.add_argument("--input", required=True, help="Input FASTQ(.gz) containing candidate reads")
    parser.add_argument("--output", required=True, help="Output FASTQ(.gz) containing rescued reads")
    parser.add_argument("--stats", required=True, help="Output JSON stats path")
    parser.add_argument("--threshold", type=float, default=0.6, help="Confidence threshold")
    parser.add_argument("--max-reads", type=int, default=5000, help="Max reads to score")
    parser.add_argument("--url", default="https://health.api.nvidia.com/v1/biology/arc/evo2-7b/generate", help="Evo2 API URL")
    parser.add_argument("--log", default=None, help="Optional log file path")
    args = parser.parse_args()

    return {
        "unclassified_fastq": args.input,
        "output_fastq": args.output,
        "output_stats": args.stats,
        "threshold": args.threshold,
        "max_reads_to_process": args.max_reads,
        "evo2_url": args.url,
        "log_path": args.log,
    }

async def get_confidence_score(session, sequence, read_id, api_url, retries=3):
    """
    Call Evo2 API to get confidence score for a sequence with retry logic.
    """
    if not NVIDIA_API_KEY:
        return {"id": read_id, "confidence": 0.0, "error": "No API Key"}

    headers = {
        "Authorization": f"Bearer {NVIDIA_API_KEY}",
        "Content-Type": "application/json"
    }

    # Use first 500bp
    analysis_seq = sequence[:500] if len(sequence) > 500 else sequence
    
    payload = {
        "sequence": analysis_seq,
        "num_tokens": 1,
        "temperature": 0.1,
        "top_k": 1,
        "enable_sampled_probs": True
    }

    for attempt in range(retries):
        try:
            async with session.post(api_url, json=payload, headers=headers, timeout=30) as response:
                if response.status == 200:
                    data = await response.json()
                    sampled_probs = data.get("sampled_probs", [])
                    if sampled_probs:
                        avg_prob = sum(sampled_probs) / len(sampled_probs)
                        return {"id": read_id, "confidence": avg_prob}
                    else:
                        return {"id": read_id, "confidence": 0.0, "error": "No probs"}
                elif response.status == 429:
                    # Rate limited - wait and retry
                    wait_time = (attempt + 1) * 2
                    await asyncio.sleep(wait_time)
                    continue
                else:
                    text = await response.text()
                    return {"id": read_id, "confidence": 0.0, "error": f"Status {response.status}: {text[:100]}"}
        except Exception as e:
            if attempt == retries - 1:
                return {"id": read_id, "confidence": 0.0, "error": str(e)}
            await asyncio.sleep(1)
            
    return {"id": read_id, "confidence": 0.0, "error": "Max retries exceeded"}

async def process_reads(unclassified_fastq, output_fastq, threshold, max_reads_to_process, api_url, log):
    log.write(f"Starting EVO2 Rescue for {unclassified_fastq}\n")
    log.write(f"Confidence threshold: {threshold}\n")
    log.flush()

    # Read candidate reads
    reads = []
    handle = _open_text_maybe_gzip(unclassified_fastq, "rt")
    
    for i, record in enumerate(SeqIO.parse(handle, "fastq")):
        if i >= max_reads_to_process:
            log.write(f"Reached limit of {max_reads_to_process} reads. Sampling...\n")
            break
        reads.append((str(record.seq).upper(), record.id))
    handle.close()

    if not reads:
        log.write("No reads to process.\n")
        return [], 0

    log.write(f"Processing {len(reads)} reads via EVO2 API...\n")
    log.flush()

    rescued_ids = set()
    confidences = []

    async with aiohttp.ClientSession() as session:
        # Process in smaller chunks to avoid hitting rate limits
        chunk_size = 10
        for i in range(0, len(reads), chunk_size):
            chunk = reads[i:i+chunk_size]
            tasks = [get_confidence_score(session, seq, rid, api_url) for seq, rid in chunk]
            results = await asyncio.gather(*tasks)
            
            for res in results:
                if "error" in res:
                    # If it's still a 429 after retries, just log it
                    log.write(f"Issue for {res['id']}: {res['error']}\n")
                
                conf = res.get("confidence", 0.0)
                confidences.append(conf)
                if conf >= threshold:
                    rescued_ids.add(res['id'])
            
            log.write(f"Processed {min(i+chunk_size, len(reads))}/{len(reads)} reads... Rescued: {len(rescued_ids)}\n")
            log.flush()
            
            # Small sleep between chunks to stay under rate limit
            await asyncio.sleep(0.5)

    # Write rescued reads
    input_handle = _open_text_maybe_gzip(unclassified_fastq, "rt")
    output_handle = _open_text_maybe_gzip(output_fastq, "wt")
    
    written = 0
    for record in SeqIO.parse(input_handle, "fastq"):
        if record.id in rescued_ids:
            SeqIO.write(record, output_handle, "fastq")
            written += 1
    
    output_handle.close()
    input_handle.close()
    
    return confidences, written

def main():
    try:
        io = _get_io_from_snakemake_or_cli()
        unclassified_fastq = io["unclassified_fastq"]
        output_fastq = io["output_fastq"]
        output_stats = io["output_stats"]
        threshold = io["threshold"]
        max_reads_to_process = io["max_reads_to_process"]
        log_path = io["log_path"]

        log = open(log_path, "w") if log_path else sys.stderr

        loop = asyncio.get_event_loop()
        confidences, written = loop.run_until_complete(
            process_reads(
                unclassified_fastq=unclassified_fastq,
                output_fastq=output_fastq,
                threshold=threshold,
                max_reads_to_process=max_reads_to_process,
                api_url=io["evo2_url"],
                log=log,
            )
        )
        
        log.write(f"Successfully rescued {written} reads.\n")
        
        # Stats
        stats = {
            "total_processed": len(confidences),
            "rescued_count": written,
            "threshold": threshold,
            "avg_confidence": sum(confidences)/len(confidences) if confidences else 0
        }
        with open(output_stats, 'w') as f:
            json.dump(stats, f, indent=2)
            
    except Exception as e:
        try:
            log.write(f"CRITICAL ERROR: {str(e)}\n")
        except Exception:
            pass
        import traceback
        try:
            log.write(traceback.format_exc())
        except Exception:
            pass
        sys.exit(1)
    finally:
        try:
            if log is not sys.stderr:
                log.close()
        except Exception:
            pass

if __name__ == "__main__":
    main()
