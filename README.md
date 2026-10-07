# Vibrion Sentinel

**Open genomic-epidemiological intelligence for *Vibrio cholerae* in Haiti.**

Haiti's 2010 cholera outbreak killed roughly 10,000 people. After elimination was declared in 2019, cholera came back in 2022. Vibrion Sentinel exists so the next resurgence is caught in public sequence data as early as possible: it watches every new *V. cholerae* genome linked to Haiti, characterizes it automatically, and interprets it against the full 2010–2022 historical record.

This is not a case-count dashboard. It is a **sequence-first early-warning system** — new public genomes in, population context out.

---

## How it works

```
NCBI SRA / ENA ── every 6 hours ──► Watcher (scripts/watch_sra.py)
                                          │ new accessions
                                          ▼
                     Lite pipeline, one GitHub Actions job per sample
                       ENA download → QC → map to 2010 Haiti reference
                       → SNP distance → toxin/virulence loci
                       → consensus genome → report.json
                                          │
              ┌───────────────────────────┼───────────────────────────┐
              ▼                           ▼                           ▼
     Supabase `sentinel_runs`      R2 archive                  GitHub artifacts
     (one JSON row / sample)   raw reads + consensus          per-run report.json
```

### 1. Watcher

`.github/workflows/sentinel-watch.yml` (schedule `17 */6 * * *`) runs `scripts/watch_sra.py`, which queries NCBI E-utilities for *Vibrio cholerae* sequencing runs in two priority tiers:

- **haiti** — any record linked to Haiti/Caribbean, any date, not yet seen
- **global** — new records from the last 7 days

New accessions are diffed against `state/seen_accessions.txt` and dispatched as a processing matrix (up to 10 per cycle, 2 concurrent jobs). Processed accessions are committed back to the state file, so nothing is ever processed twice. Transient NCBI rate-limits (HTTP 429) are retried with exponential backoff honoring the `Retry-After` header, and an optional `NCBI_API_KEY` repo secret raises the E-utilities limit from 3 to 10 req/s.

### 2. Lite pipeline (`workflow/sentinel_lite/`)

Each accession gets a fast, reference-based characterization:

- **QC** — mean depth, consensus called %, *V. cholerae* species purity (Kraken2)
- **Mapping** against 2010EL-1786, the 2010 Haiti outbreak reference (7PET lineage)
- **SNP distance** vs the reference — the primary divergence signal
- **Surveillance loci** — cholera toxin cluster (*ctxA*, *ctxB*, *zot*, *ace*), *tcpA*, *toxR*, O1 antigen genes (*wbeT*, *rfbV*), and more
- **Consensus genome** — archived for downstream phylogeny

Output is a single `report.json` per sample: machine-readable, pushed to Supabase, with a per-run artifact on the workflow run. Samples that fail QC are reported as failed, not silently dropped. Haiti-tier raw reads go to `r2://raw/haiti/{ACC}/`; QC-pass consensus genomes to `r2://consensus/{ACC}.fasta.gz`.

**Platforms and tiers (Phase 0 Nanopore branch).** The lite pipeline now runs on both Illumina and Oxford Nanopore reads, selected per-sample by `platform` in the run config (`illumina` default; `nanopore` requires `basecaller_model: fast|hac|sup` — the run refuses to start without it). Two tiers:

| Tier | What it does | When to use it |
|------|--------------|----------------|
| **lite** (Tier A) | Mapping-based: NanoPlot QC → chopper filter (Q≥9, length≥1000) → minimap2 `map-ont` → the same Kraken2 / mapping / SNP / loci / consensus path as Illumina (bwa/fastp kept for Illumina) | Every sample, every platform — the fast first pass |
| **assembly** (Tier B) | Flye (`--nano-hq`/`--nano-raw` by basecaller) → Medaka polish (basecaller-matched model) → assembly QC gate (3.8–4.4 Mb, ≤50 contigs, N50 ≥200 kb) → MLST (PubMLST *V. cholerae*) → vibecheck lineage → AMRFinderPlus (gene hits only) | Nanopore only, when lineage/AMR calls are needed |

The assembly QC gate is load-bearing: a fragmented assembly invalidates the whole sample — lineage and AMR calls are suppressed, and the failure is recorded in the report rather than silently passed through. AMR output is **gene detected/not detected only, never a susceptibility prediction** (every report carries that disclaimer verbatim). Nanopore QC thresholds are tuned for the platform (site depth 15, mean depth 30, called 85%). Every `report.json` records `platform`, `basecaller_model`, tool+database versions, and carries pipeline version `0.2.0`.

> **Decision record — Snakemake, not Nextflow (Phase 0).** The Nanopore branch extends the existing Snakemake lite pipeline; there is no migration to Nextflow in Phase 0. Rationale: the pipeline, its tests, and its CI are Snakemake-shaped; a rewrite would re-litigate every validated behavior for no Phase 0 gain. The Nextflow question reopens only if Phase 1 job orchestration forces it.

### 3. Historical cohorts (`cohorts/`)

| Cohort | Runs | Span | Source |
|--------|------|------|--------|
| `haiti-baseline.tsv` | 432 | 2010–2022 | 11 BioProjects; ENA query `tax_eq(666) AND country="Haiti*" AND instrument_platform="ILLUMINA" AND library_strategy="WGS"` |
| `backtest-2022.tsv` | 25 | Oct–Nov 2022 | PRJNA900623, the published 2022 resurgence set |

The baseline is the empirical anchor: every new genome is interpreted against what Haitian *V. cholerae* actually looked like across the outbreak (2010–11), the endemic years (2013–18), and the 2022 comeback. Construction notes and caveats: [docs/HAITI_COHORTS.md](docs/HAITI_COHORTS.md).

### 4. Back-test — would we have caught 2022 early?

`cohorts/backtest-2022.expect.yaml` encodes the published findings ([JCM 2023](https://journals.asm.org/doi/full/10.1128/jcm.00142-23); [EID 2023](https://wwwnc.cdc.gov/eid/article/29/10/23-0554_article)) as machine-checkable expectations: the 2022 isolates were toxigenic O1, homogeneous, and closely related to 2012–2019 Haitian strains. `scripts/backtest_report.py` tests the pipeline's actual outputs against those expectations. **A failure here is a finding about the pipeline and must be reported as such, not silently retuned until it passes.**

The question this project answers: *given the 2010–2018 anchor, at what point in October 2022 would the current pipeline have raised the alarm?* Note the 2022 signal is not a huge SNP distance — it is a new toxigenic O1 cluster appearing after a multi-year gap.

### 5. Phylogeny (`.github/workflows/cohort-phylo.yml`)

Core-SNP phylogeny over cohort consensus genomes — for placing new samples in the population tree rather than judging them by reference distance alone.

---

## Quick start

### Check a watcher run against the published record

```bash
# download the per-sample artifacts from a sentinel-watch run, then:
python3 scripts/backtest_report.py \
  --expect cohorts/backtest-2022.expect.yaml \
  --results-dir <artifacts dir>
```

### Force-process accessions

Actions → **Vibrion sentinel watch** → Run workflow → `force_accessions: SRR22265444,SRR22265443` (comma-separated), `force_priority: haiti`.

### Query the data

- **Supabase** `sentinel_runs`: one row per sample — accession, priority, QC status/reasons, depth, breadth, `snps_vs_7pet`, surveillance loci, and the full `report.json`.
- **R2**: raw Haiti-tier reads and QC-pass consensus genomes (see paths above).
- **GitHub**: per-run `sentinel-{ACC}` artifacts with `report.json`.

### Run the watcher search locally

```bash
python3 scripts/watch_sra.py --state state/seen_accessions.txt \
  --days 7 --max-new 10 --out-matrix /tmp/matrix.json
```

---

## The v2.0 laboratory pipeline

`workflow/Snakefile` is the original 45-rule Snakemake pipeline: deep per-sample characterization — serotype, *ctxB* allele, CTXφ integration, SXT element assembly, AMR (targeted + RGI), phenotype prediction, Pilon-polished consensus, MAFFT/FastTree phylogeny. It remains the deep-dive engine for samples the lite pipeline flags. Region configs live alongside it (`workflow/haiti_2026_config.yaml`, `workflow/haiti_resurgence_paired.yaml`, …).

```bash
# 1. Clone
git clone https://github.com/intelogroup/vibrion-sentinel.git
cd vibrion-sentinel

# 2. Create environment
conda env create -f environment.yml
conda activate vibrion

# 3. Download databases (~10GB total)
bash scripts/setup_databases.sh

# 4. Run (laboratory mode)
bash scripts/run_pipeline.sh --config workflow/test_config.yaml --cores 8
```

Results → `data/pipeline_output/MY_SAMPLE/08_comprehensive_report/surveillance_report.md`

System requirements: 8 GB RAM minimum (16 GB recommended), 2+ CPU cores, 15 GB free disk; the pipeline auto-selects a memory tier (FULL / BALANCED / BUNKER / EMERGENCY). HyenaDNA local triage setup: [docs/HYENADNA_SETUP.md](docs/HYENADNA_SETUP.md).

> AI/ML placement is deliberate: rigorous genomics and epidemiology come first. Optional model-based scoring only ever ranks already-detected candidates — it never gates, suppresses, or creates findings.

---

## Phase 1: upload ingestion (`service/`)

Field labs push raw FASTQ to Sentinel instead of waiting for public
archives. tus 1.0.0 resumable uploads in, pipeline `report.json` out.

```
lab ── tus (POST /uploads, PATCH chunks, HEAD resume) ──► ingest-api
                                                            │ job row (Supabase ingest_jobs)
                                                            ▼ R2 intake/{org}/{job}/{file}
                                                      ingest-worker ──► snakemake (sentinel_lite)
                                                            │ report.json → job row + sentinel_runs mirror
lab ◄── GET /jobs/{id} (status + report) ──────────────────┘
```

### Deploy

```bash
# 1. Create the tables (once, Supabase SQL editor):
#    sql/ingest_jobs.sql   sql/org_api_keys.sql
# 2. Create an API key for the lab (prints once — store it safely):
python -m service.make_key --org-id mirebalais --name mirebalais-lab-uploader
# 3. Export env and bring up the service:
export SUPABASE_URL=... SUPABASE_SERVICE_ROLE_KEY=... \
       R2_ENDPOINT_URL=... R2_ACCESS_KEY_ID=... R2_SECRET_ACCESS_KEY=... \
       R2_BUCKET=... PIPELINE_REFDIR=/app/data/references
docker compose up --build ingest-api ingest-worker
```

### Environment variables

| Variable | Required | Default | Purpose |
|---|---|---|---|
| `SUPABASE_URL` | yes | — | PostgREST endpoint (same project as `sentinel_runs`) |
| `SUPABASE_SERVICE_ROLE_KEY` | yes | — | service-role key; API + worker bypass RLS with it |
| `R2_ENDPOINT_URL` | yes | — | S3-compatible endpoint |
| `R2_ACCESS_KEY_ID` / `R2_SECRET_ACCESS_KEY` | yes | — | R2 credentials (env only, never baked into the image) |
| `R2_BUCKET` | yes | — | intake bucket |
| `PIPELINE_REFDIR` | yes (worker) | `/app/data/references` | reference bundle for the pipeline |
| `SNAKEMAKE_CORES` | no | `4` | worker parallelism |
| `MAX_UPLOAD_BYTES` | no | `10000000000` (10 GB) | per-file cap (multiplexed runs are 2–8 GB) |
| `WORKER_POLL_INTERVAL` | no | `10` | seconds between claim attempts |
| `REPO_DIR` | no | `/app` | repo root inside the image (locates the Snakefile) |
| `API_PORT` | no | `8000` | host port for the API |

### Upload protocol (tus 1.0.0 core)

`Upload-Metadata` keys: `filename`, `platform` (`illumina`|`nanopore`),
`basecaller_model` (`fast`|`hac`|`sup`), `tier` (`lite`|`assembly`, default
`lite`), optional `checksum` (sha256 hex, verified at completion).
Fail-closed exactly like the pipeline: `platform=nanopore` without
`basecaller_model` → 400 at creation; a parity test
(`service/tests/test_tus.py::test_validation_parity_with_pipeline`) asserts
the API and `report_lib.validate_platform_config` accept/reject the same
inputs. Auth is per-org bearer keys; one org's resources read as 404 to
another (no existence oracle). Job statuses move forward only
(`uploading → queued → running → done|failed`, `cancelled` from
`uploading`/`queued`, `failed → queued` only via explicit retry) —
enforced in code *and* by a database trigger (`sql/ingest_jobs.sql`).

### What Phase 1 does NOT build (handoffs)

- **Phase 2:** per-org encryption keys, data-deletion workflows, roles/teams/SSO,
  full RBAC. The `org_id` scoping and `org_api_keys` table are the seam it builds on.
- **Phase 3:** French/Creole PDF rendering, cloud batch/autoscaling workers
  (Phase 1 runs one local worker via `FOR UPDATE SKIP LOCKED` claims, safe to
  scale later), billing/metering.
- Production hardening not yet done: rate limiting, request logging/audit trail,
  R2 multipart streaming (Phase 1 stages chunks to local disk and PUTs once at
  completion — see `service/app/storage.py`), TLS termination.

## Roadmap

- **Strain database** — curated Postgres/PostgREST store: every Haitian/Caribbean isolate with collection date, department, clinical vs environmental origin, lineage, AMR, toxin, QC flags
- **Population analytics** — cumulative core-SNP alignment, time-scaled tree, lineage assignment, AMR/toxin trend analysis
- **Anomaly detection** — stats-first flags: novel lineage emergence, AMR acquisition, SNP-distance outliers, geographic jumps
- **Outputs** — public dashboard + monthly *"What changed in the Haitian V. cholerae population?"* brief
- **Manuscript** — methods paper + data note on the curated strain set

---

## Citation

If you use Vibrion Sentinel in published work:

> Vibrion Sentinel — open genomic-epidemiological intelligence for *Vibrio cholerae* in Haiti.
> https://github.com/intelogroup/vibrion-sentinel

Key published references for the 2022 resurgence data:
- J Clin Microbiol 2023 — https://journals.asm.org/doi/full/10.1128/jcm.00142-23
- Emerg Infect Dis 2023 — https://wwwnc.cdc.gov/eid/article/29/10/23-0554_article

---

## License

MIT License. See `LICENSE` for details.
