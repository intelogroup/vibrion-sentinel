import subprocess
import json

def get_run_from_id(sra_id):
    cmd_sum = f"curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id={sra_id}&retmode=json'"
    try:
        sum_res = subprocess.run(cmd_sum, shell=True, capture_output=True, text=True)
        sum_data = json.loads(sum_res.stdout)
        result = sum_data.get("result", {}).get(sra_id, {})
        runs = result.get("runs", "")
        print(f"ID {sra_id} Runs: {runs}")
    except:
        pass

print(" checking ID 113254...")
get_run_from_id("113254")

print("Searching CIRS 101...")
# Try adding [Organism] or searching text
subprocess.run("curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=CIRS_101&retmode=json'", shell=True)

