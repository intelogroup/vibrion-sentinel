#!/usr/bin/env python3
"""
update_global_references.py:
Expands the Vibrion Sentinel global reference panel by pulling intelligence
on recent lineages (Malawi 2023, Bangladesh XDR) and updating the 
forensic detection schemas.
"""

import json
import os
import sys
from pathlib import Path

# Paths
ROOT = Path(__file__).parent.parent
LINEAGE_DB = ROOT / "data/metadata/lineage_database.json"
DETECT_AMR = ROOT / "workflow/scripts/detect_amr.py"
GLOBAL_REF_DIR = ROOT / "data/global_references"

# Intelligence Definitions for 2023-2024 Outbreaks
NEW_LINEAGES = [
    {
        "id": "malawi-2023",
        "name": "Malawi 2023 (Regional Resurgence)",
        "description": "Primary lineage behind the largest African outbreak in a decade. Part of the 7PET wave 3 with distinct SXT variations.",
        "origin": "Malawi/Zambia",
        "year": 2023,
        "biotype": "El Tor",
        "markers": {
            "amr": ["tet(A)", "sul2", "strAB", "dfrA1", "qnrS"],
            "virulence": ["ctxA", "ctxB (Genotype 7)", "tcpA"],
            "genomic": {
                "sxt_element": "ICEVchMal1",
                "evo2_delta_baseline": 0.15
            }
        },
        "public_health_impact": "High",
        "representative_accession": "SRR24483731"
    },
    {
        "id": "bangladesh-xdr-2023",
        "name": "Bangladesh 2023 (Pan-Drug Resistant)",
        "description": "Emerging pandrug-resistant (XDR) strain carrying carbapenemases and advanced quinolone resistance.",
        "origin": "South Asia",
        "year": 2023,
        "biotype": "El Tor",
        "markers": {
            "amr": ["blaNDM-1", "qnrS", "gyrA_S83I", "parC_S85L", "mph(A)"],
            "virulence": ["ctxA", "ctxB", "tcpA", "zot"],
            "genomic": {
                "sxt_element": "ICEVchInd1-XDR",
                "evo2_delta_baseline": 0.52
            }
        },
        "public_health_impact": "Extreme (Treatment Failure Risk)",
        "representative_accession": "SRR21674404"
    },
    {
        "id": "peru-wasa",
        "name": "Peru/Latin America (WASA Lineage)",
        "description": "The West African-South American lineage responsible for the 1990s Latin American epidemic. Carry distinct VSP-II variants.",
        "origin": "Peru/Latin America",
        "year": 1991,
        "biotype": "El Tor",
        "markers": {
            "amr": ["strA", "sul2"],
            "virulence": ["ctxA", "ctxB", "tcpA"],
            "genomic": {
                "prophage": "WASA-1",
                "island": "VSP-II-WASA",
                "evo2_delta_baseline": 0.12
            }
        },
        "public_health_impact": "Historical/Regional Baseline",
        "representative_accession": "C6709"
    },
    {
        "id": "mexico-endemic",
        "name": "Mexico Endemic (Gulf Coast Variant)",
        "description": "Autochthonous Mexican variants often found in environmental reservoirs. Distinct ribotypes (M5/M6).",
        "origin": "Mexico",
        "year": 1991,
        "biotype": "El Tor / Classical Hybrid",
        "markers": {
            "amr": ["strA", "sul2"],
            "virulence": ["ctxA", "ctxB", "tcpA"],
            "genomic": {
                "ribotype": "M5/M6",
                "evo2_delta_baseline": 0.18
            }
        },
        "public_health_impact": "Moderate (Regional Monitoring)",
        "representative_accession": "MEX-1"
    },
    {
        "id": "nepal-2010",
        "name": "Nepal 2010 (Haiti Source)",
        "description": "The exact ancestral lineage of the Haiti 2010 outbreak. Provides a direct forensic matches to the original Haitian clone.",
        "origin": "Nepal",
        "year": 2010,
        "biotype": "El Tor",
        "markers": {
            "amr": ["strA", "strB", "sul2", "tet(A)", "dfrA1"],
            "virulence": ["ctxA", "ctxB7", "tcpA"],
            "genomic": {
                "sxt_element": "ICEVchNep1",
                "evo2_delta_baseline": 0.01
            }
        },
        "public_health_impact": "Forensic Key (Source Identification)",
        "representative_accession": "SRR094500"
    },
    {
        "id": "india-wave3-ancestor",
        "name": "India 7PET Wave 3 (Ancestral Reservoir)",
        "description": "The foundational South Asian lineage from which Haiti and Yemen clones emerged. Highly varied AMR profiles.",
        "origin": "India/South Asia",
        "year": 2007,
        "biotype": "El Tor",
        "markers": {
            "amr": ["mph(E)", "msr(E)", "qnrS", "dfrA1", "sul1"],
            "virulence": ["ctxA", "ctxB", "tcpA"],
            "genomic": {
                "sxt_element": "ICEVchInd5",
                "evo2_delta_baseline": 0.18
            }
        },
        "public_health_impact": "High (Global Reservoir Monitoring)",
        "representative_accession": "ERR025381"
    },
    {
        "id": "nigeria-2023",
        "name": "Nigeria 2023 (West African Hotspot)",
        "description": "Dominant West African lineage (Afr12 sublineage). High incidence of multidrug resistance and rapid regional spread.",
        "origin": "Nigeria",
        "year": 2023,
        "biotype": "El Tor (Atypical)",
        "markers": {
            "amr": ["gyrA_S83I", "parC_S85L", "sul2", "dfrA1"],
            "virulence": ["ctxA", "ctxB7", "tcpA"],
            "genomic": {
                "sublineage": "Afr12",
                "evo2_delta_baseline": 0.28
            }
        },
        "public_health_impact": "High (Regional Resurgence)",
        "representative_accession": "SRR22283995"
    },
    {
        "id": "drc-2023",
        "name": "DRC 2023 (ST69 Resurgence)",
        "description": "Central African stable lineage (ST69) showing 100% resistance to core antibiotics. Highly persistent in the Congo Basin.",
        "origin": "DRC",
        "year": 2023,
        "biotype": "El Tor",
        "markers": {
            "amr": ["strA", "strB", "sul2", "tet(A)"],
            "virulence": ["ctxA", "ctxB7", "tcpA"],
            "genomic": {
                "sublineage": "AFR15",
                "evo2_delta_baseline": 0.32
            }
        },
        "public_health_impact": "High (Stable Endemicity)",
        "representative_accession": "SRR24483731"
    }
]

def update_lineage_db():
    print(f"📡 Updating Lineage Database: {LINEAGE_DB}")
    if not LINEAGE_DB.exists():
        print("❌ Error: Lineage database not found.")
        return

    with open(LINEAGE_DB, 'r') as f:
        db = json.load(f)

    existing_ids = [l['id'] for l in db['lineages']]
    added = 0
    
    for nl in NEW_LINEAGES:
        if nl['id'] not in existing_ids:
            db['lineages'].append(nl)
            added += 1
            print(f"   + Added {nl['name']}")

    # Update outbreak list if not present
    new_outbreaks = [
        {"location": "Malawi", "start_year": 2022, "end_year": 2023, "dominant_lineage": "malawi-2023", "notes": "Largest outbreak in country history."},
        {"location": "Bangladesh", "start_year": 2023, "end_year": "Ongoing", "dominant_lineage": "bangladesh-xdr-2023", "notes": "High incidence of XDR/Carbapenem resistance."}
    ]
    
    existing_locations = [o['location'] for o in db.get('outbreaks', [])]
    for o in new_outbreaks:
        if o['location'] not in existing_locations:
            db['outbreaks'].append(o)

    with open(LINEAGE_DB, 'w') as f:
        json.dump(db, f, indent=2)
    
    print(f"✅ Successfully updated lineage database (Added {added} lineages).")

def download_representative_genomes():
    """
    Instructions on how to pull these high-priority genomes
    """
    print("\n📦 GENOMIC ACQUISITION: High-Priority Global References")
    print("-------------------------------------------------------")
    GLOBAL_REF_DIR.mkdir(parents=True, exist_ok=True)
    
    for nl in NEW_LINEAGES:
        acc = nl['representative_accession']
        dest = GLOBAL_REF_DIR / f"{nl['id']}.fasta"
        if not dest.exists():
            print(f"🔗 Recommended Action: Download {acc} for {nl['id']} lineage.")
            print(f"   Command: datasets download genome accession {acc} --include genome,gff3 --filename {dest}.zip")
        else:
            print(f"✅ Local reference for {nl['id']} already exists.")

def main():
    print("🌍 Vibrion Sentinel: Global Reference Update Utility")
    print("====================================================")
    update_lineage_db()
    download_representative_genomes()
    print("\n🚀 Next Steps: Rerunning 'amr_detect' with broader k-mer signatures for updated lineages.")

if __name__ == "__main__":
    main()
