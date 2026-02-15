import os
import json

def main():
    mapping_file = "data/serogroup_reference/serogroup_mapping.json"
    with open(mapping_file, "r") as f:
        mapping = json.load(f)
    
    # Unique serogroups
    serogroups = sorted(list(set(v["serogroup"] for v in mapping.values())))
    
    # 1. Prepare names.dmp
    # Format: taxid | name | unique name | name class |
    names = []
    names.append("1\t|\tall\t|\t\t|\tscientific name\t|")
    names.append("2\t|\tBacteria\t|\t\t|\tscientific name\t|")
    names.append("666\t|\tVibrio cholerae\t|\t\t|\tscientific name\t|")
    
    # 2. Prepare nodes.dmp
    # Format: taxid | parent taxid | rank | ...
    nodes = []
    nodes.append("1\t|\t1\t|\tno rank\t|\t\t|")
    nodes.append("2\t|\t1\t|\tsuperkingdom\t|\t\t|")
    nodes.append("666\t|\t2\t|\tspecies\t|\t\t|")
    
    serogroup_to_taxid = {}
    next_id = 9000001
    
    for sg in serogroups:
        taxid = next_id
        serogroup_to_taxid[sg] = taxid
        
        names.append(f"{taxid}\t|\t{sg}\t|\t\t|\tscientific name\t|")
        nodes.append(f"{taxid}\t|\t666\t|\tserogroup\t|\t\t|")
        
        next_id += 1
    
    taxonomy_dir = "data/kraken2_serogroup/taxonomy"
    if not os.path.exists(taxonomy_dir):
        os.makedirs(taxonomy_dir)
        
    with open(os.path.join(taxonomy_dir, "names.dmp"), "w") as f:
        f.write("\n".join(names) + "\n")
        
    with open(os.path.join(taxonomy_dir, "nodes.dmp"), "w") as f:
        f.write("\n".join(nodes) + "\n")
        
    with open("data/serogroup_reference/serogroup_to_taxid.json", "w") as f:
        json.dump(serogroup_to_taxid, f, indent=4)
        
    print(f"Created taxonomy with {len(serogroups)} serogroups.")

if __name__ == "__main__":
    main()
