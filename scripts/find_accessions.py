import subprocess
import json
import time

def search_sra(term):
    print(f"Searching SRA for: {term}")
    # esearch -db sra -query "term" | efetch -format runinfo
    # We will simulate this using curl to the API or just use the python requests if available? 
    # Actually, let's use the provided `run_command` with curl to the eutils API.
    # Base URL: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=...&retmode=json
    
    cmd_search = f"curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={term}&retmode=json'"
    
    try:
        search_res = subprocess.run(cmd_search, shell=True, capture_output=True, text=True)
        data = json.loads(search_res.stdout)
        id_list = data.get("esearchresult", {}).get("idlist", [])
        
        print(f"Found IDs: {id_list}")
        
        for sra_id in id_list[:3]: # Check first 3
            # Get run info
            # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=ID&retmode=json
            cmd_sum = f"curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id={sra_id}&retmode=json'"
            sum_res = subprocess.run(cmd_sum, shell=True, capture_output=True, text=True)
            sum_data = json.loads(sum_res.stdout)
            
            result = sum_data.get("result", {}).get(sra_id, {})
            # print(json.dumps(result, indent=2))
            
            # Extract run accession
            runs = result.get("runs", "")
            # The runs field in JSON summary is often a JSON string or structure, need to parse
            # Actually, standard esummary JSON for SRA is messy. 
            # Let's just output the text summary to grep
            print(f"Summary for {sra_id}: {runs}")
            
    except Exception as e:
        print(f"Error: {e}")

# Candidate 3
search_sra("2012EL-1410")

# Candidate 4 (Dirty Bomb)
search_sra("hc-17a1")
search_sra("hc-77a1")

# Candidate 5 (Altered El Tor)
search_sra("Vibrio cholerae Matlab variant")
search_sra("CIRS101")

# Candidate 6 (False Flag)
# Check Walters 2023 / PRJNA900623
search_sra("PRJNA900623")
