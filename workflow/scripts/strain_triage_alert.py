#!/usr/bin/env python3
import json
import os
import sys
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Strain Triage Alert System")
    parser.add_argument("--serogroup", required=True, help="Path to serogroup report JSON")
    parser.add_argument("--output", required=True, help="Path to save triage decision JSON")
    parser.add_argument("--sample", required=True, help="Sample ID")
    args = parser.parse_args()

    print(f"🕵️  Triage system analyzing sample: {args.sample}")
    
    if not os.path.exists(args.serogroup):
        print(f"❌ Error: Serogroup report not found at {args.serogroup}")
        sys.exit(1)

    with open(args.serogroup, 'r') as f:
        report = json.load(f)

def analyze_triage_decision(report, sample_id):
    serogroup = report.get("serogroup", report.get("primary_serogroup", "Unknown"))
    
    # Try multiple keys for serotype
    serotype = report.get("serotype", "Unknown")
    if serotype == "Unknown":
        serotype_details = report.get("serotype_details", {})
        if isinstance(serotype_details, dict):
            serotype = serotype_details.get("serotype_status", "Unknown")

    lineage = report.get("lineage_context", "Unknown")
    toxin = report.get("virulence_markers", {}).get("status", "Unknown")
    toxin_type = report.get("virulence_markers", {}).get("type", "NONE")

    # Define the "Safe/Standard" lineage (Haiti O1)
    # Both Ogawa and Inaba are now standard for Haiti 2022
    EXPECTED_SEROGROUP = "O1"
    
    is_atypical = False
    alert_level = "INFO"
    routing_mode = "STRICT_SURVEILLANCE"
    forensic_directives = []
    message = f"Sample confirmed as {serogroup} ({serotype}). Matches standard surveillance profile."

    # Logic Flow
    if serogroup == "O1":
        # Standard O1 (Ogawa or Inaba)
        if serotype == "Inaba (Likely)":
             message += " Note: Inaba serotype detected (wbeT mutation confirmed)."
        elif serotype == "Ogawa (Provisional)":
             message += " Note: Ogawa serotype detected."
        else:
             # O1 with unknown serotype (Rough?)
             alert_level = "WARNING"
             message = f"⚠️  WARNING: O1 Serogroup confirmed but Serotype undefined ({serotype}). Check for O-antigen loss."
        
        # 2029 Future-Proof: Check for Hybrid Virulence even in O1
        if toxin_type == "NOVC_PATHOGEN" or "chxA" in str(report):
             is_atypical = True
             alert_level = "CRITICAL"
             routing_mode = "FORENSIC_DISCOVERY"
             message = f"🚨 ALERT: DETECTED O1 HYBRID STRAIN. Contains NOVC Virulence factors (Cholix/T3SS). Pathogen evolution detected."
             forensic_directives.append("Scan for Horizontal Gene Transfer (HGT) events")
             forensic_directives.append("Verify T3SS islet integration site")
             
    elif "O139" in serogroup:
        is_atypical = True
        alert_level = "CRITICAL"
        routing_mode = "FORENSIC_DISCOVERY"
        message = f"🚨 ALERT: DETECTED O139 STRAIN ({serogroup}). Switching to FORENSIC_DISCOVERY mode."
        forensic_directives.append("Compare against Bengal 1993 archetype")
        forensic_directives.append("Scan for O139-specific wbeT/wbf clusters")
        
    elif "NON_O1" in serogroup or "Non-O1" in serogroup:
         # Non-O1 Logic
         if toxin == "VIRULENT":
            is_atypical = True
            if toxin_type == "NOVC_PATHOGEN":
                alert_level = "CRITICAL"
                message = f"🚨 CRITICAL: DETECTED PATHOGENIC Non-O1 strain (T3SS/Cholix). Switching to FORENSIC_DISCOVERY mode."
                forensic_directives.append("Scan for T3SS island integrity")
            else:
                alert_level = "WARNING"
                message = f"⚠️  WARNING: DETECTED VIRULENT Non-O1/Non-O139 strain. Switching to FORENSIC_DISCOVERY mode."
                forensic_directives.append("Scan for hybrid virulence elements")
            routing_mode = "FORENSIC_DISCOVERY"
         else:
             # Benign Non-O1
             message = "Non-O1/Non-O139 strain detected (Benign/Environmental)."
             
    elif lineage != "Unknown" and "Haiti" not in lineage:
        # Fallback to lineage check if available
        is_atypical = True
        alert_level = "ELEVATED"
        routing_mode = "FORENSIC_DISCOVERY"
        message = f"⚠️  ELEVATED: Divergent lineage detected: {lineage}. Switching to FORENSIC_DISCOVERY mode."
        forensic_directives.append("Perform full-spectrum delta-anomaly against multiple archetypes")

    decision = {
        "sample_id": sample_id,
        "is_atypical": is_atypical,
        "alert_level": alert_level,
        "routing_mode": routing_mode,
        "forensic_directives": forensic_directives,
        "message": message,
        "serogroup": serogroup,
        "lineage": lineage,
        "proceed_to_discovery": is_atypical,
        "timestamp": report.get("timestamp")
    }
    return decision

def main():
    parser = argparse.ArgumentParser(description="Strain Triage Alert System")
    parser.add_argument("--serogroup", required=True, help="Path to serogroup report JSON")
    parser.add_argument("--output", required=True, help="Path to save triage decision JSON")
    parser.add_argument("--sample", required=True, help="Sample ID")
    args = parser.parse_args()

    print(f"🕵️  Triage system analyzing sample: {args.sample}")
    
    if not os.path.exists(args.serogroup):
        print(f"❌ Error: Serogroup report not found at {args.serogroup}")
        sys.exit(1)

    with open(args.serogroup, 'r') as f:
        report = json.load(f)
        
    decision = analyze_triage_decision(report, args.sample)
    
    # Print summary to stdout
    print("-" * 60)
    if decision["is_atypical"]:
        print(f"[{decision['alert_level']}] {decision['message']}")
        print(f"💡 Routing to: {decision['routing_mode']}")
        for directive in decision['forensic_directives']:
            print(f"   - {directive}")
    else:
        print(f"✅ Sample verified: {decision['lineage']}")
    print("-" * 60)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(decision, f, indent=2)

    print(f"Triage decision saved to {args.output}")

if __name__ == "__main__":
    main()
