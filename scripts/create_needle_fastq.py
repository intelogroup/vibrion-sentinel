import os
import random
from Bio import SeqIO
import gzip

# Paths
VIBRION_DIR = "/Users/kalinovdameus/Developer/Vibrion"
NEEDLE_FASTA = os.path.join(VIBRION_DIR, "data/references/2010EL-1786.fasta")
HAYSTACK_DIR = os.path.join(VIBRION_DIR, "data/kraken2_library")

# Background genomes
HAYSTACK_FILES = [
    os.path.join(HAYSTACK_DIR, "other_vibrio/parahaemolyticus.fasta"),
    os.path.join(HAYSTACK_DIR, "other_vibrio/vulnificus.fasta"),
    os.path.join(HAYSTACK_DIR, "novc_environmental/archaea.fasta"),
    os.path.join(HAYSTACK_DIR, "novc_environmental/fungi.fasta")
]

OUTPUT_FASTQ = os.path.join(VIBRION_DIR, "data/validation/needle_haystack_test.fastq.gz")
os.makedirs(os.path.dirname(OUTPUT_FASTQ), exist_ok=True)

READ_LEN = 150
TOTAL_READS = 100000
NEEDLE_RATIO = 0.01  # 1%
NEEDLE_COUNT = int(TOTAL_READS * NEEDLE_RATIO)
HAYSTACK_COUNT = TOTAL_READS - NEEDLE_COUNT

def get_random_read(seq, read_len, record_id, count):
    reads = []
    seq_len = len(seq)
    for i in range(count):
        start = random.randint(0, seq_len - read_len - 1)
        subseq = seq[start:start+read_len]
        # Simulate base quality (all high quality for now)
        qual = "I" * read_len
        read_id = f"@{record_id}_{i+start}"
        reads.append(f"{read_id}\n{subseq}\n+\n{qual}\n")
    return reads

def main():
    print(f"Generating synthetic 'Needle-in-a-Haystack' test...")
    print(f"Target: 1% V. cholerae 2010EL-1786 ({NEEDLE_COUNT} reads)")
    print(f"Background: 99% environmental microbial DNA ({HAYSTACK_COUNT} reads)")

    all_reads = []

    # 1. Generate Needle reads
    needle_record = list(SeqIO.parse(NEEDLE_FASTA, "fasta"))[0]
    all_reads.extend(get_random_read(needle_record.seq, READ_LEN, "VC_2010", NEEDLE_COUNT))
    print(f"✓ Generated {NEEDLE_COUNT} needle reads")

    # 2. Generate Haystack reads
    reads_per_background = HAYSTACK_COUNT // len(HAYSTACK_FILES)
    for bg_file in HAYSTACK_FILES:
        try:
            bg_records = list(SeqIO.parse(bg_file, "fasta"))
            # Use the biggest chromosome/contig
            bg_record = sorted(bg_records, key=lambda x: len(x.seq), reverse=True)[0]
            bg_name = os.path.basename(bg_file).split(".")[0]
            all_reads.extend(get_random_read(bg_record.seq, READ_LEN, f"BG_{bg_name}", reads_per_background))
            print(f"✓ Generated {reads_per_background} reads from {bg_name}")
        except Exception as e:
            print(f"Error reading {bg_file}: {e}")

    # Shuffle reads
    random.shuffle(all_reads)

    # Write to gzipped FASTQ
    with gzip.open(OUTPUT_FASTQ, "wt") as f:
        f.writelines(all_reads)

    print(f"\n✅ Synthetic sample created: {OUTPUT_FASTQ}")
    print(f"Total reads: {len(all_reads)}")

if __name__ == "__main__":
    main()
