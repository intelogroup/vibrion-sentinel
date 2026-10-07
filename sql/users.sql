-- Vibrion Sentinel: human user accounts (Phase 3 portal).
--
-- Password auth (simpler ops than magic-link email delivery; documented
-- tradeoff in the README). Passwords are argon2 hashes — never MD5/SHA.
-- API keys (sql/org_api_keys.sql) keep working unchanged for
-- machine-to-machine access.
--
-- Run this once in the Supabase SQL editor (or `psql`).

create table if not exists users (
    id            uuid primary key default gen_random_uuid(),
    email         text not null unique,      -- lowercased at write time
    password_hash text not null,             -- argon2id hash
    created_at    timestamptz not null default now()
);

create index if not exists users_email_idx on users (email);

alter table users enable row level security;

drop policy if exists "users read access" on users;
create policy "users read access"
    on users for select
    using (true);
