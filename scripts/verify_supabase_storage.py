import requests
import os
import sys

SUPABASE_URL = "https://ndrfttboavxgtsvatcim.supabase.co"
SUPABASE_SERVICE_KEY = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6Im5kcmZ0dGJvYXZ4Z3RzdmF0Y2ltIiwicm9sZSI6InNlcnZpY2Vfcm9sZSIsImlhdCI6MTc2ODk2MzE3NCwiZXhwIjoyMDg0NTM5MTc0fQ.jf8QIskOlDwthKAXhdHe2Rf1Z7AbYha8tKGd_UUo1mM"

headers = {
    "Authorization": f"Bearer {SUPABASE_SERVICE_KEY}",
    "apikey": SUPABASE_SERVICE_KEY
}

def list_buckets():
    response = requests.get(f"{SUPABASE_URL}/storage/v1/bucket", headers=headers)
    if response.status_code == 200:
        buckets = response.json()
        print(f"Buckets found: {[b['name'] for b in buckets]}")
        return buckets
    else:
        print(f"Failed to list buckets: {response.status_code} {response.text}")
        return []

def list_files(bucket_name, path=""):
    response = requests.post(
        f"{SUPABASE_URL}/storage/v1/object/list/{bucket_name}",
        headers=headers,
        json={"prefix": path, "limit": 100, "offset": 0, "sortBy": {"column": "name", "order": "asc"}}
    )
    if response.status_code == 200:
        files = response.json()
        print(f"\nFiles in bucket '{bucket_name}' (path: '{path}'):")
        for f in files:
            metadata = f.get('metadata') or {}
            size = metadata.get('size', 'unknown')
            print(f"  - {f['name']} ({size} bytes)")
    else:
        print(f"Failed to list files in {bucket_name}: {response.status_code} {response.text}")

def check_tables():
    # PostgREST API
    tables = ["samples", "genomes", "variants", "analysis_results"]
    print("\nChecking PostgREST Tables:")
    for table in tables:
        response = requests.get(
            f"{SUPABASE_URL}/rest/v1/{table}?select=count",
            headers={**headers, "Prefer": "count=exact"}
        )
        if response.status_code == 200:
            count = response.headers.get("Content-Range", "0-0/0").split("/")[-1]
            print(f"  - {table}: {count} rows")
        elif response.status_code == 404:
            print(f"  - {table}: Not found")
        else:
            print(f"  - {table}: Error {response.status_code}")

def main():
    buckets = list_buckets()
    for bucket in buckets:
        list_files(bucket['name'])
    
    check_tables()

if __name__ == "__main__":
    main()
