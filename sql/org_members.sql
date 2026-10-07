-- Vibrion Sentinel: human membership in orgs (Phase 3 portal).
--
-- Roles reuse the API-key role set: admin | member | viewer. A human's
-- effective permission in an org is their membership role. API keys are
-- unaffected (machine-to-machine keeps working).
--
-- Run this once in the Supabase SQL editor (or `psql`).

create table if not exists org_members (
    user_id     uuid not null references users (id) on delete cascade,
    org_id      text not null,
    role        text not null default 'member'
                check (role in ('admin', 'member', 'viewer')),
    created_at  timestamptz not null default now(),
    primary key (user_id, org_id)
);

create index if not exists org_members_org_id_idx on org_members (org_id);

alter table org_members enable row level security;

drop policy if exists "org_members read access" on org_members;
create policy "org_members read access"
    on org_members for select
    using (true);
