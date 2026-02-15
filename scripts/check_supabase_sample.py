import asyncio
import sys
from pathlib import Path

# Add project root to path
sys.path.append(str(Path(__file__).resolve().parents[1]))

from backend.services.supabase_storage import SupabaseStorageService

async def check_sample(prefix):
    print(f"🔍 Checking Supabase for: {prefix}")
    try:
        storage = SupabaseStorageService()
        files = await storage.list_genomes(prefix=prefix)
        
        if not files:
            print(f"❌ No files found matching '{prefix}' in Supabase.")
        else:
            print(f"✅ Found {len(files)} files:")
            for f in files:
                metadata = f.get('metadata') or {}
                size = metadata.get('size', 0)
                print(f"  - {f.get('name')} ({size / (1024*1024):.2f} MB)")
                
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    prefix = sys.argv[1] if len(sys.argv) > 1 else ""
    asyncio.run(check_sample(prefix))
