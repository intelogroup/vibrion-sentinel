-- Vibrion Sentinel: per-run metrics history.
--
-- One row per processed SRA/ENA accession, upserted by the sentinel-watch
-- GitHub Actions workflow after each successful pipeline run. Raw FASTQ is
-- NOT stored here (see R2 archive for Haiti-tier raw reads) -- this table
-- is the queryable summary: coverage, SNP counts, surveillance-loci calls.
--
-- Run this once in the Supabase SQL editor (or `psql` against the project's
-- connection string) to create the table.

create table if not exists sentinel_runs (
    id                  bigint generated always as identity primary key,
    accession           text not null unique,
    priority            text,                    -- 'haiti' | 'global' | null (forced/manual)
    run_at              timestamptz not null default now(),

    -- workflow/sentinel_lite/Snakefile rule `report`
    pipeline_version    text,
    reference            text,

    reads_raw           integer,
    reads_qc_passed     integer,
    reads_vibrio_kept   integer,

    kraken_db           text,
    v_cholerae_reads    integer,
    v_cholerae_fraction numeric,

    mapped_reads        integer,
    mean_depth          numeric,
    breadth_pct         numeric,

    snps_vs_7pet         integer,
    consensus_length     integer,
    consensus_called_pct numeric,   -- % of reference not masked (depth >= min_site_depth)

    qc_status            text,      -- 'pass' | 'fail': filter on this before using variants/loci
    qc_reasons           jsonb,

    surveillance_loci    jsonb,   -- {locus: {call: present/partial/absent, present, mean_depth, breadth_pct}}
    report                jsonb,   -- full report.json, kept verbatim for anything not modeled above

    created_at           timestamptz not null default now()
);

-- for tables created before the QC columns existed
alter table sentinel_runs add column if not exists consensus_called_pct numeric;
alter table sentinel_runs add column if not exists qc_status text;
alter table sentinel_runs add column if not exists qc_reasons jsonb;

create index if not exists sentinel_runs_qc_status_idx on sentinel_runs (qc_status);
create index if not exists sentinel_runs_priority_idx on sentinel_runs (priority);
create index if not exists sentinel_runs_run_at_idx on sentinel_runs (run_at desc);

-- Row Level Security: table is written only by the pipeline's service-role
-- key (bypasses RLS by default). Enable + add a read policy if you want
-- e.g. an anon/authenticated dashboard to query this directly.
alter table sentinel_runs enable row level security;

drop policy if exists "sentinel_runs read access" on sentinel_runs;
create policy "sentinel_runs read access"
    on sentinel_runs for select
    using (true);
