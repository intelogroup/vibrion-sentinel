#!/usr/bin/env python3
import json
import argparse
import subprocess
import os
import tempfile
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path

# --- Constants & Coordinates (2010EL-1786 / CP003069.1) ---
# Used for VCF Heterozygosity Checks
COORDS = {
    "ctxB": {"chrom": "CP003069.1", "start": 1041237, "end": 1041612}, # 1-based rough range
    "rstR": {"chrom": "CP003069.1", "start": 1047218, "end": 1047556},
    "tcpA": {"chrom": "CP003069.1", "start": 367950, "end": 368624},
    # "wbeT" is handled in track_wbeT_mutations.py, but we could add it here for completeness
}

def get_ctxb_genotype(fasta_path):
    """
    Genotype 1 (Classical): 39-His (H), 68-Thr (T)
    Genotype 3 (El Tor): 39-Tyr (Y), 68-Ile (I)
    """
    # VC1456 ctxB gene coordinates in 2010EL-1786 (Chromosome 1)
    # We will use targeted BLAST to extract the gene
    query_gene = "VC1456" 
    
    # Representative ctxB gene sequence for extraction (Classical reference)
    ctxb_ref = "ATGATTAAGATTATTTGCGTCGCTCAGACGGGATTTGTTAGGCTGGATGGCTCAAAATATTGCATACCAATTCTATTTCTACAGACAGAGTGAGTACTTTGACAGGTTTGAGGGCAATTACACAGATGCGGAGCGAGAAAGACGAACGTATGGGTTAAA"
    
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(">ctxB_ref\n" + ctxb_ref + "\n")
        temp_query = f.name
        
    # Use -task blastn-short with relaxed word size for fragmented genome components
    cmd = [
        "blastn", 
        "-query", temp_query, 
        "-subject", fasta_path, 
        "-outfmt", "6 sseq sstrand", 
        "-task", "blastn-short",
        "-evalue", "1e-5"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if os.path.exists(temp_query): os.unlink(temp_query)

    is_classical_39 = False
    is_eltor_39 = False
    aa_39, aa_68 = "?", "?"

    if result.stdout.strip():
        # Extract the first hit (highest score)
        lines = result.stdout.strip().split("\n")
        hit_parts = lines[0].split("\t")
        ctxb_seq = hit_parts[0].replace("-", "").upper()
        
        # Analyze Nucleotide markers directly if translation is risky
        is_classical_39 = "AATCACACA" in ctxb_seq 
        is_eltor_39 = "AATTACACA" in ctxb_seq
        
        # Try translation
        try:
            if len(ctxb_seq) >= 150:
                protein = str(Seq(ctxb_seq).translate())
                aa_39 = protein[38] if len(protein) > 38 else "?"
                aa_68 = protein[67] if len(protein) > 67 else "?"
            else:
                if is_classical_39: aa_39 = "H"
                elif is_eltor_39: aa_39 = "Y"
        except Exception:
            pass
    
    # FALLBACK: Direct search in the FASTA file if BLAST was inconclusive or residues missing
    if aa_39 == "?" and not is_classical_39 and not is_eltor_39:
        # Check for Classical H39 signature
        cmd_c = ["grep", "-c", "AATCACACA", fasta_path]
        res_c = subprocess.run(cmd_c, capture_output=True, text=True).stdout.strip()
        if res_c.isdigit() and int(res_c) > 0:
            is_classical_39 = True
            aa_39 = "H"
        else:
            # Check reverse complement
            cmd_cr = ["grep", "-c", "TGTGTGATT", fasta_path]
            res_cr = subprocess.run(cmd_cr, capture_output=True, text=True).stdout.strip()
            if res_cr.isdigit() and int(res_cr) > 0:
                is_classical_39 = True
                aa_39 = "H"
        
        # Check for El Tor Y39 signature
        if not is_classical_39:
            cmd_e = ["grep", "-c", "AATTACACA", fasta_path]
            res_e = subprocess.run(cmd_e, capture_output=True, text=True).stdout.strip()
            if res_e.isdigit() and int(res_e) > 0:
                is_eltor_39 = True
                aa_39 = "Y"
            else:
                cmd_er = ["grep", "-c", "TGTGTAATT", fasta_path]
                res_er = subprocess.run(cmd_er, capture_output=True, text=True).stdout.strip()
                if res_er.isdigit() and int(res_er) > 0:
                    is_eltor_39 = True
                    aa_39 = "Y"

    genotype = "Unknown"
    
    # Simple assignment for voting
    allele = "Unknown"
    if is_classical_39:
        genotype = "Genotype 1 (Classical)"
        allele = "Classical"
    elif is_eltor_39:
        genotype = "Genotype 3 (El Tor)"
        allele = "El Tor"
    
    return {
        "genotype": genotype,
        "allele": allele,
        "aa_39": aa_39,
        "aa_68": aa_68,
        "ctxb_gene_present": is_classical_39 or is_eltor_39 or aa_39 != "?"
    }

# Helper for sequence searching (handles RC)
def search_in_genome(genome_seq, queries):
    """
    Search for a list of query sequences in the genome (and its RC).
    Returns True if ANY query is found.
    """
    seq_upper = genome_seq.upper()
    rc_seq = str(Seq(seq_upper).reverse_complement())
    
    for q in queries:
        q_upper = q.upper()
        if q_upper in seq_upper or q_upper in rc_seq:
            return True
    return False

def get_rstr_allele(fasta_path):
    """
    rstR^Classical vs rstR^ElTor detection.
    Different alleles of the rstR repressor gene in the CTX phage.
    """
    # Load genome once (lazy approach, could be optimized to pass seq)
    try:
        with open(fasta_path, 'r') as f:
            # Simple fasta parsing: ignore headers > check sequences
            # Better: concat all sequence lines
            lines = [l.strip() for l in f if not l.startswith(">")]
            genome_seq = "".join(lines)
    except:
        return {"allele": "Error reading file", "status": "ERROR"}

    # Specific markers
    # Validated from 2010EL-1786 (Haitian El Tor)
    markers = {
        "rstR_Classical": ["ATGATAAACGATGCGCGTTGGATTACGACGATG"], 
        "rstR_ElTor": ["ATGAAGATAAAAGAAAGGCTAGCCAACCAA"] # Validated 2010EL-1786 Start
    }
    
    detected = []
    
    for allele, seqs in markers.items():
        if search_in_genome(genome_seq, seqs):
            detected.append(allele)
                
    if "rstR_Classical" in detected and "rstR_ElTor" in detected:
         return {"allele": "See details", "status": "MIXED", "detected": detected}
    elif "rstR_Classical" in detected:
         return {"allele": "Classical", "status": "CLASSICAL_PHAGE_REPRESSOR"}
    elif "rstR_ElTor" in detected:
         return {"allele": "El Tor", "status": "EL_TOR_PHAGE_REPRESSOR"}
    else:
         return {"allele": "Unknown/Missing", "status": "ABSENT"}

def get_tcpA_allele(fasta_path):
    """
    tcpA Classical vs El Tor discrimination.
    """
    try:
        with open(fasta_path, 'r') as f:
            lines = [l.strip() for l in f if not l.startswith(">")]
            genome_seq = "".join(lines)
    except:
        return {"allele": "Error", "status": "ERROR"}

    markers = {
        "tcpA_Classical": ["GACTTTTGAGAT"], # Classical specific
        "tcpA_ElTor": ["CGTGCGATTGATTCGCAGAAT"] # Validated 2010EL-1786 Unique Region
    }
    
    detected = []
    
    for allele, seqs in markers.items():
        if search_in_genome(genome_seq, seqs):
            detected.append(allele)
            
    if "tcpA_Classical" in detected:
        return {"allele": "Classical", "status": "CLASSICAL_PILI"}
    elif "tcpA_ElTor" in detected:
        return {"allele": "El Tor", "status": "EL_TOR_PILI"}
    else:
        return {"allele": "Unknown/Missing", "status": "ABSENT"}

def get_hylA_allele(fasta_path):
    """
    hylA (Hemolysin) Indel Check.
    """
    try:
        with open(fasta_path, 'r') as f:
            lines = [l.strip() for l in f if not l.startswith(">")]
            genome_seq = "".join(lines)
    except:
        return {"allele": "Error", "status": "ERROR"}
        
    # 11bp: TTTTTAGCATT (Example deleted region in Classical)
    eltor_11bp_insert = "TTTTTAGCATT" 
    
    if search_in_genome(genome_seq, [eltor_11bp_insert]):
        return {"allele": "El Tor (Functional)", "status": "HEMOLYSIN_ACTIVE"}
    
    return {"allele": "Classical/Unknown (Deletion or Absent)", "status": "HEMOLYSIN_INACTIVE"}

def check_7pet_markers(fasta_path):
    """
    Check for 7th Pandemic (El Tor) specific islands: VSP-I, VSP-II, RS1 (rstC).
    """
    try:
        with open(fasta_path, 'r') as f:
            lines = [l.strip() for l in f if not l.startswith(">")]
            genome_seq = "".join(lines)
    except:
        return {"RS1_rstC": "ERROR", "VSP_I": "ERROR", "VSP_II": "ERROR"}

    markers = {
        "RS1_rstC": "ATGAACAAATCTCAAGAAATGGCTATCAAT", 
        "VSP_I": "TTAGCGTTGGTCGAGCGCATA", 
        "VSP_II": "ATGACACAAACAAATCAAACT" 
    }
    
    found = {}
    for name, seq in markers.items():
         if search_in_genome(genome_seq, [seq]):
             found[name] = "PRESENT"
         else:
             found[name] = "ABSENT"
            
    return found

def check_region_heterozygosity(vcf_path, region_name):
    """
    Scans VCF for heterozygous calls in the specified gene region.
    Returns details if mixed.
    """
    if not vcf_path or not os.path.exists(vcf_path):
        return None

    target = COORDS.get(region_name)
    if not target: return None
    
    chrom = target["chrom"]
    start = target["start"]
    end = target["end"]
    
    mixed_sites = []
    
    try:
        if vcf_path.endswith(".gz"):
            opener = gzip.open(vcf_path, 'rt')
        else:
            opener = open(vcf_path, 'r')
            
        with opener as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.split('\t')
                if len(parts) < 10: continue
                
                v_chrom = parts[0]
                v_pos = int(parts[1])
                
                if v_chrom != chrom: continue
                if v_pos < start: continue
                if v_pos > end: break # Sorted VCF assumption
                
                # Check AF
                info = parts[7]
                fmt = parts[8]
                smp = parts[9]
                
                af = None
                if "AF=" in info:
                    for tag in info.split(';'):
                        if tag.startswith("AF="):
                            af = float(tag.split('=')[1])
                            break
                elif "AD" in fmt:
                    try:
                        ad_idx = fmt.split(':').index("AD")
                        ad_str = smp.split(':')[ad_idx]
                        counts = [int(x) for x in ad_str.split(',')]
                        if sum(counts) > 0:
                            af = counts[1] / sum(counts)
                    except: pass
                
                if af is not None and 0.15 <= af <= 0.85: # 15-85% is definitely mixed
                    mixed_sites.append(f"Pos {v_pos}: AF={af:.2f}")
                    
    except Exception:
        pass

    if mixed_sites:
        return {"status": "MIXED", "sites": mixed_sites, "count": len(mixed_sites)}
    return {"status": "PURE"}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Fasta file")
    parser.add_argument("--vcf", help="VCF file for heterozygosity check")
    parser.add_argument("--output", help="Output JSON")
    args = parser.parse_args()
    
    # 1. Run Individual Marker Checks
    ctxb = get_ctxb_genotype(args.input)
    rstr = get_rstr_allele(args.input)
    tcpA = get_tcpA_allele(args.input)
    hylA = get_hylA_allele(args.input)
    arch_7pet = check_7pet_markers(args.input)
    
    # 2. VCF Heterozygosity Checks
    mixtures = {}
    if args.vcf:
        for gene in ["ctxB", "rstR", "tcpA"]:
            res = check_region_heterozygosity(args.vcf, gene)
            if res and res["status"] == "MIXED":
                mixtures[gene] = res
    
    # 3. Vote & Confidence
    score_classical = 0.0
    score_eltor = 0.0
    
    # Core Markers (Heavy Weight)
    if "Classical" in ctxb["allele"]: score_classical += 1
    elif "El Tor" in ctxb["allele"]: score_eltor += 1
    
    if "Classical" in rstr["allele"]: score_classical += 1
    elif "El Tor" in rstr["allele"]: score_eltor += 1
    
    if "Classical" in tcpA["allele"]: score_classical += 1
    elif "El Tor" in tcpA["allele"]: score_eltor += 1

    # Supportive Markers (Lower Weight)
    # hylA downgraded to 0.3
    if "Classical" in hylA["allele"]: score_classical += 0.3
    elif "El Tor" in hylA["allele"]: score_eltor += 0.3
    
    # 7PET Boosters
    # Each present 7PET marker adds 0.2 to El Tor score
    for m, status in arch_7pet.items():
        if status == "PRESENT":
            score_eltor += 0.2

    # 4. Biotype Determination
    biotype_backbone = "Unknown"
    
    # House markers (rstR + tcpA) define the backbone
    backbone_cl = 0
    backbone_et = 0
    
    if "Classical" in rstr["allele"]: backbone_cl += 1
    elif "El Tor" in rstr["allele"]: backbone_et += 1
    
    if "Classical" in tcpA["allele"]: backbone_cl += 1
    elif "El Tor" in tcpA["allele"]: backbone_et += 1
    
    if backbone_cl > backbone_et:
        biotype_backbone = "Classical"
    elif backbone_et > backbone_cl:
        biotype_backbone = "El Tor"
    
    # 5. Discordance Matrix
    # Compare ctxB, rstR, tcpA
    markers_list = [ctxb["allele"], rstr["allele"], tcpA["allele"]]
    unique_calls = set(m for m in markers_list if m not in ["Unknown/Missing", "Unknown"])
    
    is_discordant = False
    discordance_type = "None"
    
    if "Classical" in unique_calls and "El Tor" in unique_calls:
        is_discordant = True
        discordance_type = "Mixed Backbone/Toxin"
    
    # 6. Final Variant Classification
    variant_pathotype = "Unknown"
    
    if biotype_backbone == "El Tor":
        if "Classical" in ctxb["genotype"]:
             variant_pathotype = "Altered El Tor (Mozambique-like)"
             if is_discordant: discordance_type = "Altered El Tor Pattern"
        else:
             variant_pathotype = "Typical El Tor (7th Pandemic)"
    elif biotype_backbone == "Classical":
         variant_pathotype = "Classical (6th Pandemic)"
         
    # Handle Hikojima / Mixtures explicitly
    if mixtures:
        variant_pathotype += " / MIXED (Hikojima?)"
        is_discordant = True
        discordance_type = "Heterozygous Mixture"

    # 7. Confidence Score (0-100)
    # Simple logic: Ratio of matching markers to total checked
    total_markers = 3.0 # ctxB, rstR, tcpA
    matches = 0.0
    if biotype_backbone == "El Tor":
        matches = score_eltor # approximate
    elif biotype_backbone == "Classical":
        matches = score_classical
        
    # Cap at 100
    confidence = min(100, int((matches / 3.0) * 100))
    if is_discordant and not "Altered" in variant_pathotype:
        confidence -= 20 # Penalty for unexplained discordance
    
    # Altered El Tor is a accepted discordance, so high confidence is possible
    if "Altered" in variant_pathotype:
        confidence = 95 # High confidence in the Altered state

    report = {
        "biotype": biotype_backbone,
        "variant": variant_pathotype,
        "confidence_score": confidence,
        "discordance": {
            "is_discordant": is_discordant,
            "type": discordance_type,
            "matrix": {
                "ctxB": ctxb["allele"],
                "rstR": rstr["allele"],
                "tcpA": tcpA["allele"]
            }
        },
        "mixtures": mixtures,
        "markers": {
            "ctxb": ctxb,
            "rstr": rstr,
            "tcpA": tcpA,
            "hylA": hylA,
            "7pet_arch": arch_7pet
        },
        "scores": {"Classical": score_classical, "El Tor": score_eltor}
    }
    
    if args.output:
        with open(args.output, "w") as f:
            json.dump(report, f, indent=4)
    else:
        print(json.dumps(report, indent=4))

if __name__ == "__main__":
    main()
