-- Migration 004 (Phase 3): allow running -> queued for bounded infra retries.
--
-- The worker retries INFRA errors only (download/R2 failures), max 3, with
-- backoff. The bound lives in code (service/worker.py); the database permits
-- the transition so the retry path is legal. Pipeline errors never rewind:
-- they go running -> failed, and only the explicit POST /jobs/{id}/retry
-- moves failed -> queued.
--
-- Legal after this migration:
--   uploading -> queued|failed|cancelled
--   queued    -> running|failed|cancelled
--   running   -> done|failed|queued   (queued = infra retry only, max 3)
--   failed    -> queued              (explicit retry only)

alter table ingest_jobs
    add column if not exists retry_count integer not null default 0;

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
        when 'running'   then array['done', 'failed', 'queued']
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
