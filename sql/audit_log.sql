-- Vibrion Sentinel: append-only audit log (Phase 2).
--
-- Logged: key.created, key.revoked, dek.rotated, job.deleted, org.purged.
-- Append-only is enforced three ways:
--   1. The service code exposes only INSERT (db.audit) and SELECT
--      (db.list_audit) — there is no update/delete method to call.
--   2. The trigger below rejects UPDATE and DELETE at the database level.
--   3. This file documents the intent: the audit trail is evidence, and
--      evidence must not be editable through the service.
--
-- (The service_role key bypasses RLS/policies; the trigger is the real
-- backstop. A database superuser could still drop the trigger — that is
-- outside the service's threat model and would itself be visible in the
-- database's own logs.)
--
-- Run once in the Supabase SQL editor (or psql).

create table if not exists audit_log (
    id           uuid primary key default gen_random_uuid(),
    org_id       text not null,
    actor_key_id uuid,              -- which API key acted (null only for system actions)
    action       text not null,     -- key.created | key.revoked | dek.rotated |
                                   -- job.deleted | org.purged
    target       text not null,     -- key id / job id / org id acted upon
    detail       jsonb not null default '{}',
    at           timestamptz not null default now()
);

create index if not exists audit_log_org_at_idx
    on audit_log (org_id, at desc);
create index if not exists audit_log_action_idx
    on audit_log (action);

alter table audit_log enable row level security;

drop policy if exists "audit_log read access" on audit_log;
create policy "audit_log read access"
    on audit_log for select using (true);

-- Append-only: reject UPDATE and DELETE outright.
create or replace function audit_log_no_mutation()
returns trigger as $$
begin
    raise exception 'audit_log is append-only: % not allowed', TG_OP;
    return null;
end;
$$ language plpgsql;

drop trigger if exists audit_log_no_mutation on audit_log;
create trigger audit_log_no_mutation
    before update or delete on audit_log
    for each row execute function audit_log_no_mutation();
