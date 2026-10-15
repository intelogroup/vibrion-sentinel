# HL2601 regression fixture — V. cholerae ST2040 (novel ST case)

Strain HL2601: clinical non-O1/non-O139 *V. cholerae* from a bloodstream
infection, Xiamen, China (Frontiers in Microbiology 2026). Its 7-locus MLST
profile matched nothing in PubMLST; the curator assigned novel ST2040.
Before this fixture, Sentinel's typing layer reduced that outcome to
`mlst_st: null` — indistinguishable from a failed typing run.

## Files

- `genomic.fasta.gz` — the deposited assembly, byte-for-byte the NCBI
  Datasets genome fasta for **GCF_059628145.1** (GenBank GCA_059628145.1,
  ASM5962814v1; BioProject PRJNA1479141, BioSample SAMN60942877).
- `mlst.tsv` — tseemann/mlst-shaped output for this genome: ST `-` with all
  seven allele columns called (adk 51, gyrB 41, mdh 4, metE 79, pntA 72,
  purM 1, pyrC 83). Allele numbers were verified 2026-10-15 by exact
  sequence match of each locus extracted from the assembly against the
  live PubMLST vcholerae allele database — every allele is individually
  known; the *combination* is the novel part.
- `marker_refs.fasta` — ctxA, ctxB, tcpA, rtxA, hlyA, hcp gene sequences
  from the 2010EL-1786 reference at the `loci_2010EL-1786.bed`
  coordinates, used by the regression test to derive loci calls.

## Content verification (2026-10-15, before bundle/test entry)

Label claims checked against sequence content, per repo rule:

| Check | Method | Result |
|---|---|---|
| organism | NCBI Datasets report for the accession | *Vibrio cholerae*, tax_id **666**, strain HL2601 |
| assembly length | fasta stats | **4,297,354 bp**, 47 records, N50 324,224 — in the 3.8–4.3 Mb range; passes the repo assembly QC gate |
| ctxAB | 31-mer containment of 2010EL-1786 ctxA/ctxB in the assembly; PGAP GFF | containment **0.000 / 0.000**; not annotated → **absent** |
| tcpA | same | containment **0.000**; not annotated → **absent** |
| rtxA | same | containment **0.679**; PGAP-annotated → **present** |
| hlyA | same | containment **0.821**; PGAP-annotated → **present** |
| T6SS | PGAP GFF; hcp 31-mer containment | `hcp-2` + `vasH`/vas cluster annotated → **present** (containment of the reference hcp-1 sequence is 0.227 — the HL2601 copy is a diverged paralog; exact-k-mer containment understates it, the annotation is the stronger signal) |
| MLST | per-locus exact match vs PubMLST allele DB | 7/7 loci located, all exact known alleles; combination has no ST (paper: curator assigned ST2040). Live profile-table re-check is login-gated by PubMLST for post-2024 data, so the assignment rests on the paper + the verified allele profile. |

Expected Sentinel classification of this genome: **V. cholerae,
unassigned ST**, non-toxigenic (ctxAB absent), markers tcpA absent /
rtxA present / hlyA present / T6SS present.
