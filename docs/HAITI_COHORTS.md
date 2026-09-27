# Haiti cohorts

The empirical backbone of Vibrion Sentinel: every new genome the watcher processes is interpreted against what Haitian *Vibrio cholerae* actually looked like from the 2010 outbreak through the 2022 resurgence.

## `cohorts/haiti-baseline.tsv` — the anchor

432 Illumina WGS runs from 11 BioProjects (PRJNA188139, PRJNA245985, PRJNA255987, PRJNA266293, PRJNA339415, PRJNA400505, PRJNA510624, PRJNA590944, PRJNA723557, PRJNA900623, PRJNA903489), built with the ENA portal query:

```
tax_eq(666) AND country="Haiti*" AND instrument_platform="ILLUMINA" AND library_strategy="WGS"
```

Collection-year spread:

| Years | Runs | Era |
|-------|------|-----|
| 2010–2011 | 61 | Outbreak |
| 2012 | 3 | |
| 2013–2018 | 305 | Endemic persistence |
| 2019–2020 | 0 | Elimination declared 2019 — near-absence of sequences is itself data |
| 2021 | 2 | |
| 2022 | 59 | Resurgence |

Columns: `run_accession, sample_accession, study_accession, country, collection_date, first_public, instrument_platform, library_strategy, read_count`.

### Caveats

- **Depth heterogeneity.** Read counts span ~200k to ~31M. Low-depth samples may fail QC or yield partial loci calls; depth is a property of the data, not of the strain.
- **Two undated records** (SRR23510629, SRR23510635, PRJNA266293) have no `collection_date`. They are usable for population context but must be excluded from any time-ordered anchor analysis.
- **Geography is coarse.** Most records say only "Haiti"; a subset names Port-au-Prince. Department/commune-level mapping is future curation work.
- **Clinical vs environmental** is not yet encoded in the manifest. PRJNA266293 mixes both.

## `cohorts/backtest-2022.tsv` — the test

25 runs from PRJNA900623, collected 2022-10-03 through 2022-11-21: the published 2022 resurgence set. The first two (SRR22265444, 2022-10-03; SRR22265443, 2022-10-04) are low-depth and may fail QC — the early-detection question is about the earliest **QC-pass** sample that produces a signal, not the earliest sample on the calendar.

## `cohorts/backtest-2022.expect.yaml` — the contract

Encodes the published findings ([JCM 2023](https://journals.asm.org/doi/full/10.1128/jcm.00142-23); [EID 2023](https://wwwnc.cdc.gov/eid/article/29/10/23-0554_article)) as machine-checkable expectations over `scripts/backtest_report.py` output:

- **QC pass rate ≥ 0.80** — these are published, deliberately-sequenced clinical isolates; a low pass rate indicts the pipeline, not the data.
- **Toxigenic O1 loci** (*ctxA*, *ctxB*, *zot*, *ace*, *tcpA*, *toxR*, *wbeT*, *rfbV*) present in ≥ 90% of QC-pass samples.
- **SNPs vs 7PET**: median ≤ 400, IQR ≤ 150 — "homogeneous" tested as a tight spread, not an exact value.

Explicitly **not** claimed: reproducing the published Bayesian ancestral-origin analysis. A mapping-based pipeline with no assembly cannot do that; agreement on gene content and SNP spread is supporting evidence, not reproduction.

**A failure here is a finding about the pipeline and must be reported as such, not silently retuned until it passes.**

## The anchor design

Back-test logic: hold the 2010–2018 baseline as the empirical anchor, then walk the October–November 2022 samples in collection-date order and ask at each step whether the current pipeline would have raised an alarm. The 2022 signal is *not* a huge SNP distance — published work places these isolates close to 2012–2019 Haitian strains. The detectable event is a **new toxigenic O1 cluster appearing after a multi-year near-absence of sequences**, i.e. post-gap recurrence with cluster tightness, judged against the anchor's SNP-vs-7PET envelope and the cohort phylogeny — not against the fixed reference alone.
