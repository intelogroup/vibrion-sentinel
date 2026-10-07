-- Migration 003: deleted_at tombstone on sentinel_runs (Phase 2).
--
-- When a job is deleted (DELETE /jobs/{id}), its sentinel_runs mirror row
-- is KEPT and annotated with deleted_at rather than removed: population-
-- level aggregates (lineage trends, AMR frequencies) must not silently
-- rewrite history when a source sample is deleted. Queries over "live"
-- data filter WHERE deleted_at IS NULL; historical aggregates don't.
--
-- Run once in the Supabase SQL editor (or psql). Safe to re-run.

alter table sentinel_runs
    add column if not exists deleted_at timestamptz;  -- null = live row

create index if not exists sentinel_runs_deleted_at_idx
    on sentinel_runs (deleted_at) where deleted_at is not null;
