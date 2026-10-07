-- Vibrion Sentinel: org invites (Phase 3 portal).
--
-- Admin invites someone by email; the invite token is single-use and expires
-- after 7 days. The token_hash is sha256 hex; the plaintext goes to the
-- invitee out of band (no email delivery in Phase 3 — the admin hands it
-- over). Accepting requires a logged-in account.
--
-- Run this once in the Supabase SQL editor (or `psql`).

create table if not exists org_invites (
    id          uuid primary key default gen_random_uuid(),
    org_id      text not null,
    email       text not null,               -- intended recipient (advisory)
    role        text not null default 'member'
                check (role in ('admin', 'member', 'viewer')),
    token_hash  text not null unique,        -- sha256 hex of the invite token
    created_by  text not null default 'portal',
    created_at  timestamptz not null default now(),
    expires_at  timestamptz not null,
    used_at     timestamptz                  -- null = unused
);

create index if not exists org_invites_token_hash_idx on org_invites (token_hash);
create index if not exists org_invites_org_id_idx on org_invites (org_id);

alter table org_invites enable row level security;

drop policy if exists "org_invites read access" on org_invites;
create policy "org_invites read access"
    on org_invites for select
    using (true);
