import os
import json
from Bio import SeqIO

def main():
    mapping_file = "data/serogroup_reference/serogroup_mapping.json"
    with open(mapping_file, "r") as f:
        mapping = json.load(f)
    
    taxid_file = "data/serogroup_reference/serogroup_to_taxid.json"
    with open(taxid_file, "r") as f:
        sg_to_taxid = json.load(f)
    
    output_dir = "data/kraken2_serogroup/library/added"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    print(f"Processing {len(mapping)} records...")
    for acc, info in mapping.items():
        sg = info["serogroup"]
        taxid = sg_to_taxid.get(sg)
        if not taxid:
            print(f"Warning: No taxid for serogroup {sg} (accession {acc})")
            continue
            
        src_file = info["file"]
        # Ensure path is absolute if needed, but it should work relative from root
        full_src_path = os.path.join("/Users/kalinovdameus/Developer/Vibrion", src_file)
        
        if not os.path.exists(full_src_path):
            # Try just relative if absolute fails (e.g. if root differs)
            full_src_path = src_file
            
        record = SeqIO.read(full_src_path, "fasta")
        # Kraken2 format: >accession|kraken:taxid|taxid
        record.id = f"{acc}|kraken:taxid|{taxid}"
        record.description = f"V. cholerae {sg} cluster"
        
        dest_file = os.path.join(output_dir, f"{acc}_{sg}.fasta")
        SeqIO.write(record, dest_file, "fasta")
        
    print("Library preparation complete.")

if __name__ == "__main__":
    main()
