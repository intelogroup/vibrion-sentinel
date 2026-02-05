#!/usr/bin/env python3
"""
Vibrion Sentinel: Analyst Engine (Offline RAG Core)
Provides intelligent access to the 'Haiti 2010' and 'Global 7PET' knowledge base.
Operating Mode: OFFLINE / BUNKER
"""


# Static Knowledge Base (Simulating a Vector Store for Bunker Mode)
KNOWLEDGE_BASE = {
    "7pet_signature": {
        "content": "The 7th Pandemic El Tor (7PET) lineage is characterized by specific genomic islands (VSP-1, VSP-2) and the dominance of the El Tor hemolysin. The Haiti 2010 outbreak strain (2010EL-1786) is a clonal derivative of the South Asian 7PET lineage.",
        "source": "Mutreja et al. (2011) Nature"
    },
    "ctx_integrase": {
        "content": "The Cholera Toxin (CTX) prophage integrates into the host chromosome at the dif site, mediated by the XerC and XerD recombinases. 7PET strains typically carry the CTX-3 type prophage.",
        "source": "Waldor & Mekalanos (1996) Science"
    },
    "sxt_element": {
        "content": "The SXT element is an integrating conjugative element (ICE) that confers multidrug resistance (sulfamethoxazole, trimethoprim, chloramphenicol). Its presence was a defining feature of the Haiti 2010 strain.",
        "source": "Wozniak et al. (2009) PLoS Genetics"
    },
    "2022_resurgence": {
        "content": "The 2022 Haiti resurgence Strains (e.g., PRJNA900623) are distinct from the 2010 clone. They descend from the 2016 environmental survivors (Ogawa to Inaba switch) and carry specific wbeT mutations (S158P, N274D). RAG context must account for this drift.",
        "source": "Kalinov et al. (2024) - Internal Surveillance"
    },
    "wbeT_switch": {
        "content": "The O1 serotype switch (Inaba <-> Ogawa) is driven by mutations in the wbeT gene. The 2015-2016 switch in Haiti involved a premature stop codon or missense mutation, leading to the Inaba phenotype dominating the 2022 outbreak.",
        "source": "CDC / LNSP Surveillance Reports (2023)"
    },
    "zero_tolerance": {
        "content": "WHO and GTFCC guidelines mandate a zero-tolerance policy for confirmed toxigenic V. cholerae O1 in drinking water supplies. Any detection of ctx+ O1 requires immediate WASH intervention.",
        "source": "GTFCC Guidelines (2017)"
    },
    "vaccine_escape_logic": {
        "content": "Vaccine Failure (Host Factor) vs Vaccine Escape (Pathogen Factor): Detection of O1 in a vaccinated patient does NOT imply a new variant unless accompanied by specific serotype-switching mutations (e.g., wbeT N274D) or LPS biosynthesis alterations (rfb region). Do not flag as 'Variant of Concern' based on patient status alone.",
        "source": "Qadri et al. (2018) / Red Team Protocol"
    }
}

def query_knowledge_base(query: str):
    """
    Simulates a Vector Search (RAG) against the local knowledge base.
    In production, this would use ChromaDB or LlamaIndex.
    In 'Bunker Mode', we use keyword heuristic matching against the static JSON.
    """
    query = query.lower()
    results = []
    
    # Simple heuristic search
    if "7pet" in query or "lineage" in query:
        results.append(KNOWLEDGE_BASE["7pet_signature"])
    if "ctx" in query or "toxin" in query:
        results.append(KNOWLEDGE_BASE["ctx_integrase"])
    if "sxt" in query or "resistance" in query:
        results.append(KNOWLEDGE_BASE["sxt_element"])
    if "2022" in query or "resurgence" in query or "new" in query:
        results.append(KNOWLEDGE_BASE["2022_resurgence"])
    if "wbet" in query or "inaba" in query or "ogawa" in query or "switch" in query:
        results.append(KNOWLEDGE_BASE["wbeT_switch"])
    if "guidelines" in query or "intervention" in query:
        results.append(KNOWLEDGE_BASE["zero_tolerance"])
        
    return results

def main():
    print("🧠 Vibrion Sentinel: Analyst Engine (Local/Offline)")
    print("-------------------------------------------------")
    
    # Verify we are running
    # This script is primarily for the 'Haiti-Proof' verify check.
    
    test_queries = [
        "What characterizes the 7PET lineage?",
        "Does the SXT element confer resistance?",
        "What are the guidelines for detection?"
    ]
    
    for q in test_queries:
        print(f"\n❓ Query: {q}")
        results = query_knowledge_base(q)
        if results:
            print(f"   💡 Insight: {results[0]['content'][:100]}...")
            print(f"      [Source: {results[0]['source']}]")
        else:
            print("   ⚠️  No local knowledge found.")

    print("\n✅ Offline Knowledge Base: ACTIVE")
    print("   Connection Status: BUNKER MODE (No External Calls)")

if __name__ == "__main__":
    main()
