#!/usr/bin/env python3
"""Extract surveillance loci coordinates from a GenBank record into BED format.

Used once during reference building: downloads the 2010EL-1786 GenBank record,
finds virulence / serogroup / resistance features by keyword, and writes a BED
file the lite pipeline uses for per-locus coverage calls.

Usage:
    python3 scripts/extract_loci_bed.py --genbank 2010EL-1786.gb --out loci.bed
"""

import argparse
import sys

# keyword -> locus label. Matched (case-insensitive) against /gene and /product.
LOCI_KEYWORDS = {
    "ctxA": "ctxA_toxin",
    "ctxB": "ctxB_toxin",
    "tcpA": "tcpA_pilus",
    "tcpB": "tcpB_pilus",
    "wbeT": "wbeT_serogroup",
    "wbe": "wbe_cluster",
    "rfb": "rfb_cluster",
    "wbf": "wbf_cluster_O139",
    "int": "sxt_integrase",
    "tra": "sxt_conjugation",
    "gyrA": "gyrA_qrdr",
    "parC": "parC_qrdr",
    "parE": "parE_qrdr",
    "hapR": "hapR_regulator",
    "luxO": "luxO_regulator",
    "rtxA": "rtxA_toxin",
    "hlyA": "hlyA_hemolysin",
    "toxR": "toxR_regulator",
    "toxT": "toxT_regulator",
}


def parse_genbank(path):
    """Minimal GenBank feature parser. Yields (seq_id, start0, end0, strand, gene, product)."""
    seq_id = None
    features = []
    in_features = False
    current = None

    def flush():
        nonlocal current
        if current and current.get("type") in ("gene", "CDS"):
            loc = current["loc"]
            # handle complement() and ranges; skip joins spanning contigs
            strand = "-"
            if loc.startswith("complement("):
                loc = loc[len("complement("):-1]
            else:
                strand = "+"
            if "join" in loc or "<" in loc or ">" in loc:
                current = None
                return
            try:
                start_s, end_s = loc.split("..")
                start, end = int(start_s), int(end_s)
            except ValueError:
                current = None
                return
            features.append((seq_id, start - 1, end, strand,
                             current.get("gene", ""), current.get("product", "")))
        current = None

    with open(path) as fh:
        for line in fh:
            if line.startswith("LOCUS"):
                seq_id = line.split()[1]
            elif line.startswith("FEATURES"):
                in_features = True
            elif line.startswith("ORIGIN"):
                flush()
                break
            elif in_features:
                if line[5:6] not in (" ", "") and not line.startswith(" " * 6):
                    # new feature key at column 6
                    flush()
                    parts = line[5:].split(None, 1)
                    if len(parts) == 2:
                        current = {"type": parts[0], "loc": parts[1].strip(),
                                   "gene": "", "product": ""}
                elif current is not None and line.strip().startswith("/"):
                    qual = line.strip()[1:]
                    if "=" in qual:
                        k, v = qual.split("=", 1)
                        v = v.strip('"')
                        if k == "gene":
                            current["gene"] = v
                        elif k == "product" and not current["product"]:
                            current["product"] = v
    return features


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genbank", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    features = parse_genbank(args.genbank)
    print(f"parsed {len(features)} gene/CDS features", file=sys.stderr)

    found = {}
    with open(args.out, "w") as out:
        for seq_id, start, end, strand, gene, product in features:
            hay = f"{gene} {product}".lower()
            for kw, label in LOCI_KEYWORDS.items():
                if kw.lower() in hay and label not in found:
                    found[label] = (seq_id, start, end, strand, gene or product)
                    out.write(f"{seq_id}\t{start}\t{end}\t{label}\t.\t{strand}\n")
                    break

    missing = [l for l in set(LOCI_KEYWORDS.values()) if l not in found]
    print(f"wrote {len(found)} loci to {args.out}", file=sys.stderr)
    if missing:
        print(f"WARNING: loci not found in record: {', '.join(sorted(missing))}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
