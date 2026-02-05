#!/usr/bin/env python3
"""
Shared utilities for downloading samples from ENA/SRA.
Consolidates common download operations.
"""

import subprocess
from pathlib import Path
from typing import List
import httpx


def query_ena_filereport(accession: str, timeout: float = 60.0) -> List[dict]:
    """
    Query ENA API for FASTQ download links.
    
    Args:
        accession: BioProject, sample, or run accession
        timeout: Request timeout in seconds
    
    Returns:
        List of run metadata dictionaries
    """
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = {
        "accession": accession,
        "result": "read_run",
        "format": "json",
        "fields": "run_accession,fastq_ftp,collection_date,country,scientific_name",
        "limit": 0,
    }
    
    try:
        with httpx.Client(timeout=timeout) as client:
            resp = client.get(url, params=params)
            resp.raise_for_status()
            return resp.json()
    except Exception as e:
        print(f"ENA API error for {accession}: {e}")
        return []


def extract_fastq_urls(ena_records: List[dict]) -> List[str]:
    """
    Extract FASTQ URLs from ENA records, converting FTP to HTTPS.
    
    Args:
        ena_records: List of ENA metadata dictionaries
    
    Returns:
        List of HTTPS URLs
    """
    urls = []
    for record in ena_records:
        ftp_links = (record.get("fastq_ftp") or "").split(";")
        for link in ftp_links:
            link = link.strip()
            if not link:
                continue
            # Convert FTP to HTTPS
            if link.startswith("ftp://"):
                link = link.replace("ftp://", "https://", 1)
            elif not link.startswith("http"):
                link = f"https://{link}"
            urls.append(link)
    return urls


def download_file_curl(url: str, dest_path: str, retries: int = 3) -> bool:
    """
    Download file using curl with retry logic.
    
    Args:
        url: URL to download
        dest_path: Destination file path
        retries: Number of retry attempts
    
    Returns:
        True if successful, False otherwise
    """
    dest = Path(dest_path)
    
    # Skip if already exists and has reasonable size
    if dest.exists() and dest.stat().st_size > 1024 * 1024:  # > 1 MB
        print(f"  ✓ {dest.name} exists ({dest.stat().st_size / (1024*1024):.1f} MB). Skipping.")
        return True
    
    dest.parent.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "curl",
        "-L",
        "--fail",
        "--retry", str(retries),
        "--retry-delay", "2",
        url,
        "-o", str(dest),
    ]
    
    print(f"  ⬇️  Downloading {dest.name}...")
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        print(f"  ✓ Completed {dest.name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  ✗ Download failed for {dest.name}: {e}")
        return False


def download_sample_from_ena(accession: str, output_dir: str, retries: int = 3) -> bool:
    """
    Download FASTQ files for a sample from ENA.
    
    Args:
        accession: Sample or run accession
        output_dir: Output directory for FASTQ files
        retries: Number of retry attempts per file
    
    Returns:
        True if at least one file downloaded successfully
    """
    print(f"🔍 Querying ENA for {accession}...")
    records = query_ena_filereport(accession)
    
    if not records:
        print(f"  ⚠️  No data found for {accession}")
        return False
    
    urls = extract_fastq_urls(records)
    if not urls:
        print(f"  ⚠️  No FASTQ URLs found for {accession}")
        return False
    
    print(f"  Found {len(urls)} FASTQ file(s)")
    
    success_count = 0
    for url in urls:
        filename = url.split("/")[-1]
        dest_path = Path(output_dir) / filename
        if download_file_curl(url, str(dest_path), retries):
            success_count += 1
    
    return success_count > 0
