-- Vibrion Sentinel: per-org API keys for the upload service (Phase 1).
--
-- Keys are SHA-256 hashed at rest (key_hash); the plaintext is shown ONCE at
-- creation and never stored. Auth: `Authorization: Bearer <key>`.
-- Full RBAC/roles are Phase 2; this table is deliberately minimal.
--
-- Run this once in the Supabase SQL editor (or `psql`) to create the table.
-- Create keys with:  python -m service.make_key --org-id <org> --name <label>
-- (prints the plaintext key once; store it somewhere safe).

create table if not exists org_api_keys (
    id          uuid primary key default gen_random_uuid(),
    org_id      text not null,
    key_hash    text not null unique,   -- sha256 hex of the bearer key
    name        text,                   -- human label, e.g. "mirebalais-lab-uploader"
    created_at  timestamptz not null default now(),
    revoked_at  timestamptz            -- null = active
);

create index if not exists org_api_keys_org_id_idx on org_api_keys (org_id);
create index if not exists org_api_keys_key_hash_idx on org_api_keys (key_hash);

-- Service-role writes only; enable RLS with an open read policy mirroring
-- the convention in sql/ingest_jobs.sql (service_role bypasses RLS anyway).
alter table org_api_keys enable row level security;

drop policy if exists "org_api_keys read access" on org_api_keys;
create policy "org_api_keys read access"
    on org_api_keys for select
    using (true);
