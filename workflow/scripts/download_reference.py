#!/usr/bin/env python3
"""
Download 2010EL-1786 reference genome from Supabase.
This is the Haiti baseline strain needed for minimap2 alignment.
"""
import asyncio
import sys
from pathlib import Path

# Add root to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from backend.services.supabase_storage import SupabaseStorageService


async def download_reference():
    """Download reference genome from Supabase"""
    service = SupabaseStorageService()
    
    # Full storage path in Supabase (includes bucket name)
    storage_path = "genomes/reference/2010EL-1786.fasta"
    
    # Local destination
    local_path = Path(__file__).parent.parent.parent / "data" / "references" / "2010EL-1786.fasta"
    local_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"📥 Downloading {storage_path} from Supabase...")
    print(f"   Destination: {local_path}")
    
    # Download using Supabase storage service
    content = await service.download_genome(storage_path, output_path=None)
    
    if content is None:
        print(f"❌ Failed to download {storage_path}")
        return False
    
    # Write to local file
    with open(local_path, "wb") as f:
        f.write(content)
    
    # Verify size
    size_mb = local_path.stat().st_size / (1024 * 1024)
    print(f"✅ Downloaded {size_mb:.1f} MB")
    
    # Quick validation
    with open(local_path, "r") as f:
        first_line = f.readline().strip()
        if first_line.startswith(">"):
            print(f"✅ Valid FASTA format: {first_line}")
        else:
            print("❌ Invalid FASTA format")
            return False
    
    return True


if __name__ == "__main__":
    success = asyncio.run(download_reference())
    sys.exit(0 if success else 1)
