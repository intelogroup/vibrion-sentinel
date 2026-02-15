import subprocess
import json
import re

def search_bioproject(bp):
    print(f"Searching {bp}...")
    cmd = f"curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={bp}&retmode=json&retmax=50'"
    try:
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        data = json.loads(res.stdout)
        ids = data.get("esearchresult", {}).get("idlist", [])
        
        for i in ids[:20]:
            cmd_sum = f"curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id={i}&retmode=json'"
            s_res = subprocess.run(cmd_sum, shell=True, capture_output=True, text=True)
            s_data = json.loads(s_res.stdout)
            result = s_data.get("result", {}).get(i, {})
            title = result.get("exp_xml", {}).get("Summary", {}).get("Title", str(result))
            # Just dump generic info
            print(f"ID {i}: {result.get('runs')}")
    except Exception as e:
        print(e)

# search_bioproject("PRJNA189036")
# Candidate 5 might be SRR772844 ? 

# Let's search for "Altered El Tor" specifically
# cmd = "curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=Altered+El+Tor+Bangladesh&retmode=json'"
# subprocess.run(cmd, shell=True)
