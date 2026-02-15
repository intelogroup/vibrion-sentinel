
import os

LIBRARY_PATH = "/Users/kalinovdameus/Developer/Vibrion/data/kraken2_library"
MAP_FILE = "/Users/kalinovdameus/Developer/Vibrion/data/kraken2_haiti_custom/seqid2taxid.map"

mapping = {
    "haiti_2010_lineage": 666,
    "novc_environmental": 666,
    "v_mimicus": 674,
    "other_vibrio": 670, # default to parahaemolyticus
}

def get_taxid(parent_dir, filename):
    if "parahaemolyticus" in filename: return 670
    if "vulnificus" in filename: return 672
    return mapping.get(parent_dir, 666)

with open(MAP_FILE, "w") as out:
    for root, dirs, files in os.walk(LIBRARY_PATH):
        parent = os.path.basename(root)
        for f in files:
            if f.endswith(".fasta") and os.path.getsize(os.path.join(root, f)) > 0:
                taxid = get_taxid(parent, f)
                with open(os.path.join(root, f), "r") as fasta:
                    for line in fasta:
                        if line.startswith(">"):
                            # Get the first word as the ID
                            seq_id = line[1:].split()[0]
                            out.write(f"{seq_id}\t{taxid}\n")

print(f"Created map file with local identifiers at {MAP_FILE}")
