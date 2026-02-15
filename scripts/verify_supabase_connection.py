import os
import requests
from dotenv import load_dotenv

def verify_storage():
    load_dotenv()
    
    url = os.getenv("SUPABASE_URL")
    key = os.getenv("SUPABASE_SERVICE_ROLE_KEY") # Use service role for admin access
    
    if not url or not key:
        print("❌ Error: Missing credentials in .env")
        return

    print(f"Connecting to Supabase: {url}")
    
    headers = {
        "apikey": key,
        "Authorization": f"Bearer {key}"
    }
    
    # 1. Check Storage Buckets
    try:
        response = requests.get(f"{url}/storage/v1/bucket", headers=headers)
        if response.status_code == 200:
            buckets = response.json()
            print(f"✅ Storage Connected! Found {len(buckets)} buckets.")
            for b in buckets:
                print(f"   - Bucket: {b['name']} (Public: {b['public']})")
                
            # If no buckets, create a 'genomic-data' bucket to verify write access
            if not buckets:
                print("   Creating 'genomic-data' bucket for testing...")
                create_resp = requests.post(
                    f"{url}/storage/v1/bucket", 
                    headers=headers, 
                    json={"name": "genomic-data", "public": False}
                )
                if create_resp.status_code == 200:
                     print("   ✅ Created 'genomic-data' bucket successfully.")
                else:
                     print(f"   ❌ Failed to create bucket: {create_resp.text}")
        else:
            print(f"❌ Failed to list buckets: {response.status_code} - {response.text}")
            
    except Exception as e:
        print(f"❌ Connection Error: {e}")

    # 2. Check Database (PostgREST)
    try:
        # Just check health/connectivity by listing a non-existent table or root
        response = requests.get(f"{url}/rest/v1/", headers=headers)
        # 200 OK with list of routes/tables is expected
        if response.status_code == 200:
             print("✅ Database (PostgREST) Connected.")
             tables = response.json()
             if tables:
                 print("   Tables found:")
                 for t in tables:
                      # Structure depends on root endpoint config, usually returns OpenAPI or table list
                      # If it's the root, it might return links.
                      pass
        else:
             # It might be 404 if no tables exist yet?
             print(f"ℹ️ Database response: {response.status_code}")
             
    except Exception as e:
        print(f"❌ DB Check Error: {e}")

if __name__ == "__main__":
    verify_storage()
