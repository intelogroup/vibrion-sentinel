-- Vibrion Sentinel: upload-ingestion job state (Phase 1).
--
-- One row per tus upload. The API owns transitions up to `queued`;
-- the worker owns transitions out of it. Statuses move forward only,
-- enforced BOTH in application code (service/app/jobs.py) and here by the
-- trigger below (defense in depth: a stray SQL client cannot rewind a job).
--
-- Run this once in the Supabase SQL editor (or `psql` against the project's
-- connection string) to create the table.

create table if not exists ingest_jobs (
    id              uuid primary key default gen_random_uuid(),
    org_id          text not null,
    filename        text not null,
    size_bytes      bigint,                     -- declared Upload-Length; null until known
    sha256          text,                       -- hex digest over received bytes, set at completion
    client_checksum text,                       -- optional sha256 declared by the client in Upload-Metadata
    upload_offset   bigint not null default 0,  -- bytes received so far (tus resume point)
    platform        text not null,              -- illumina | nanopore
    basecaller_model text,                      -- fast | hac | sup; required when platform=nanopore
    tier            text not null default 'lite',
    status          text not null default 'uploading',
                    -- uploading | queued | running | done | failed | cancelled
    status_reason   text,                       -- failure/cancel detail; log tail on pipeline failure
    r2_key          text,                       -- intake/{org_id}/{job_id}/{filename}
    report          jsonb,                      -- report.json verbatim when done
    created_at      timestamptz not null default now(),
    updated_at      timestamptz not null default now(),

    constraint ingest_jobs_platform_chk
        check (platform in ('illumina', 'nanopore')),
    constraint ingest_jobs_tier_chk
        check (tier in ('lite', 'assembly')),
    constraint ingest_jobs_basecaller_chk
        check (platform <> 'nanopore' or (basecaller_model is not null and basecaller_model in ('fast', 'hac', 'sup'))),
    constraint ingest_jobs_assembly_chk
        check (tier <> 'assembly' or platform = 'nanopore')
);

create index if not exists ingest_jobs_org_id_idx on ingest_jobs (org_id);
create index if not exists ingest_jobs_status_idx on ingest_jobs (status);
create index if not exists ingest_jobs_created_at_idx on ingest_jobs (created_at desc);

-- Forward-only status discipline, enforced in the database.
-- Legal: uploading -> queued|failed|cancelled
--        queued    -> running|failed|cancelled
--        running   -> done|failed
--        failed    -> queued            (explicit retry only, via POST /jobs/{id}/retry)
create or replace function enforce_ingest_job_transition()
returns trigger as $$
declare
    allowed text[];
begin
    if old.status = new.status then
        return new;
    end if;
    allowed := case old.status
        when 'uploading' then array['queued', 'failed', 'cancelled']
        when 'queued'    then array['running', 'failed', 'cancelled']
        when 'running'   then array['done', 'failed']
        when 'failed'    then array['queued']
        else array[]::text[]
    end;
    if not (new.status = any (allowed)) then
        raise exception 'illegal ingest_jobs status transition % -> %',
            old.status, new.status;
    end if;
    return new;
end;
$$ language plpgsql;

drop trigger if exists ingest_jobs_transition_guard on ingest_jobs;
create trigger ingest_jobs_transition_guard
    before update of status on ingest_jobs
    for each row execute function enforce_ingest_job_transition();

-- updated_at maintenance.
create or replace function touch_ingest_jobs_updated_at()
returns trigger as $$
begin
    new.updated_at = now();
    return new;
end;
$$ language plpgsql;

drop trigger if exists ingest_jobs_touch_updated_at on ingest_jobs;
create trigger ingest_jobs_touch_updated_at
    before update on ingest_jobs
    for each row execute function touch_ingest_jobs_updated_at();

-- Atomic claim for the worker pool. Called via PostgREST RPC
-- (POST /rest/v1/rpc/claim_next_ingest_job); the FOR UPDATE SKIP LOCKED
-- makes it safe with multiple workers. Returns the claimed row, or no rows.
create or replace function claim_next_ingest_job()
returns setof ingest_jobs as $$
    with c as (
        select id
        from ingest_jobs
        where status = 'queued'
        order by created_at
        limit 1
        for update skip locked
    )
    update ingest_jobs j
    set status = 'running', updated_at = now()
    from c
    where j.id = c.id
    returning j.*;
$$ language sql;

-- Row Level Security: written only by the service-role key (the API and the
-- worker both use it, and RLS is bypassed by service_role). Enable + add a
-- read policy if a dashboard ever needs direct reads.
alter table ingest_jobs enable row level security;

drop policy if exists "ingest_jobs read access" on ingest_jobs;
create policy "ingest_jobs read access"
    on ingest_jobs for select
    using (true);
