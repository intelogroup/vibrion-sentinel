from Bio import Entrez, SeqIO
Entrez.email = "researcher@example.com"
ids = ["LC594803", "LC594804", "LC594805", "LC594806", "LC594807"]
handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
records = list(SeqIO.parse(handle, "genbank"))
for r in records:
    for f in r.features:
        if f.type == "source":
            serovar = f.qualifiers.get("serovar", ["N/A"])[0]
            print(f"ID: {r.id}, Serovar: {serovar}")
