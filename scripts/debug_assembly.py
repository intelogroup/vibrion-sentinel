from Bio import Entrez
Entrez.email = "vibrion-surveillance@example.com"
handle = Entrez.esearch(db="assembly", term="Vibrio cholerae Haiti 2011")
record = Entrez.read(handle)
for uid in record["IdList"]:
    handle = Entrez.esummary(db="assembly", id=uid)
    summary = Entrez.read(handle)
    asm = summary["DocumentSummarySet"]["DocumentSummary"][0]
    acc = asm.get("AssemblyAccession")
    name = asm.get("AssemblyName")
    gb = asm.get("FtpPath_GenBank")
    rs = asm.get("FtpPath_RefSeq")
    print(f"UID: {uid}, Accession: {acc}, Name: {name}, GenBank: {gb}, RefSeq: {rs}")
