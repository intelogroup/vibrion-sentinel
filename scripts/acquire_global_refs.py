import os
import sys
import json
import subprocess
import argparse
from pathlib import Path

# Setup paths relative to script
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "global_references"
LINEAGE_DB = ROOT / "data" / "metadata" / "lineage_database.json"

# Ensure data directory exists
DATA_DIR.mkdir(parents=True, exist_ok=True)

class GenomicAcquirer:
    def __init__(self, dry_run=False):
        self.dry_run = dry_run
        
    def log(self, msg):
        print(f"📡 [ACQUIRE] {msg}")

    def download_sra(self, accession, target_name):
        """Fetch raw reads and perform quick assembly."""
        self.log(f"Fetching SRA reads for {accession} ({target_name})...")
        if self.dry_run:
            return
            
        target_dir = DATA_DIR / f"{target_name}_raw"
        target_dir.mkdir(exist_ok=True)
        
        # Use prefetch to pull data
        subprocess.run(["prefetch", "--max-size", "50G", accession], cwd=target_dir, check=True)
        
        # Extract FASTQ
        self.log(f"Extracting FASTQ for {accession}...")
        # Check if fastq already exists to avoid re-dumping
        if not list(target_dir.glob("*.fastq")):
            subprocess.run(["fasterq-dump", accession, "--split-files", "--threads", "4", "--outdir", str(target_dir)], cwd=target_dir, check=True)
        
        # Find FASTQ files
        fastqs = list(target_dir.glob("*.fastq"))
        if not fastqs:
            self.log(f"❌ Error: No FASTQ files found in {target_dir}")
            return
            
        r1 = target_dir / f"{accession}_1.fastq"
        r2 = target_dir / f"{accession}_2.fastq"
        single = target_dir / f"{accession}.fastq"
        
        out_dir = DATA_DIR / f"{target_name}_assembly"
        
        # Check if we should try to continue or start fresh
        spades_params = out_dir / "params.txt"
        
        if spades_params.exists():
            self.log(f"Attempting to continue existing assembly for {target_name}...")
            # Resume only requires -o
            cmd = ["spades.py", "--continue", "-o", str(out_dir)]
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError:
                self.log(f"⚠️  Resume failed for {target_name}. Restarting from scratch...")
                import shutil
                shutil.rmtree(out_dir)
                out_dir.mkdir()
                # Command for fresh run
                cmd = ["spades.py", "--careful", "-o", str(out_dir), "-t", "4", "-m", "16"]
                if r1.exists() and r2.exists():
                    cmd += ["-1", str(r1), "-2", str(r2)]
                elif single.exists():
                    cmd += ["-s", str(single)]
                else:
                    cmd += ["-s", str(fastqs[0])]
                subprocess.run(cmd, check=True)
        else:
            # Fresh run
            cmd = ["spades.py", "--careful", "-o", str(out_dir), "-t", "4", "-m", "16"]
            
            if r1.exists() and r2.exists():
                self.log(f"Detected paired-end reads for {accession}")
                cmd += ["-1", str(r1), "-2", str(r2)]
            elif single.exists():
                self.log(f"Detected single-end reads for {accession}")
                cmd += ["-s", str(single)]
            else:
                # Try to find any fastq and use as single
                self.log(f"Warning: Accession files don't match standard patterns. Using {fastqs[0].name} as single-end.")
                cmd += ["-s", str(fastqs[0])]
                
            subprocess.run(cmd, check=True)
        
        # Copy scaffold to global references
        scaffold = out_dir / "scaffolds.fasta"
        if scaffold.exists():
            final_fasta = DATA_DIR / f"{target_name}.fasta"
            subprocess.run(["cp", str(scaffold), str(final_fasta)])
            self.log(f"Successfully assembled {target_name}.fasta")
            
            # Cleanup raw reads to save space
            self.log(f"Cleaning up raw reads for {target_name}...")
            for fq in fastqs:
                fq.unlink()
            # Also remove SRA folder if exists
            sra_dir = target_dir / accession
            if sra_dir.exists():
                 import shutil
                 shutil.rmtree(sra_dir)
            self.log(f"Successfully assembled {target_name}.fasta")
        else:
            self.log(f"❌ Error: Assembly failed for {target_name}")

    def download_ftp(self, ftp_url, target_name):
        """Direct download for pre-assembled genomes."""
        self.log(f"Downloading pre-assembled genome from {ftp_url}...")
        if self.dry_run:
            return
            
        target_file = DATA_DIR / f"{target_name}.fasta"
        subprocess.run(["curl", "-o", str(target_file), ftp_url])
        self.log(f"Downloaded {target_name}.fasta")

    def run_acquisition(self, target_lineages=None):
        with open(LINEAGE_DB, 'r') as f:
            db = json.load(f)
            
        for lineage in db['lineages']:
            lineage_id = lineage['id']
            if target_lineages and lineage_id not in target_lineages:
                continue
                
            accession = lineage.get('representative_accession')
            # Check if we already have it
            if (DATA_DIR / f"{lineage_id}.fasta").exists():
                self.log(f"Reference for {lineage_id} already exists. Skipping.")
                continue
                
            if not accession or accession == "Pending" or "MOCK" in accession:
                self.log(f"No valid accession for {lineage_id}. Skipping.")
                continue

            # Logic to determine if we can use FTP or need SRA
            # For simplicity, we'll try SRA for SRR/ERR accessions
            if accession.startswith(("SRR", "ERR", "DRR")):
                self.download_sra(accession, lineage_id)
            else:
                # Placeholder for direct GenBank downloads if needed
                self.log(f"Unsupported accession type {accession} for automated download.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Vibrion Global Genomic Acquisition")
    parser.add_argument("--dry-run", action="store_true", help="Log actions without downloading")
    parser.add_argument("--lineage", nargs="+", help="Specific lineage IDs to acquire")
    args = parser.parse_args()
    
    acquirer = GenomicAcquirer(dry_run=args.dry_run)
    acquirer.run_acquisition(target_lineages=args.lineage)
