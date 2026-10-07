-- Migration 002: roles on API keys (Phase 2).
--
-- Existing keys were created before roles existed, when every key was
-- implicitly all-powerful. They become 'admin' so nothing breaks; new keys
-- default to 'member' (least privilege) and the role is chosen explicitly
-- at creation (POST /org/keys, service.make_key --role).
--
-- The backfill keys off NULL (not the 'member' value) so re-running this
-- migration never promotes keys created after it.
--
-- Run once in the Supabase SQL editor (or psql). Safe to re-run.

alter table org_api_keys
    add column if not exists role text;

-- Backfill: keys predating roles (role IS NULL) were fully privileged.
update org_api_keys set role = 'admin' where role is null;

alter table org_api_keys alter column role set default 'member';
alter table org_api_keys alter column role set not null;

-- Enforce the role vocabulary at the database level too (the API validates
-- as well; defense in depth).
do $$
begin
    if not exists (
        select 1 from pg_constraint where conname = 'org_api_keys_role_chk'
    ) then
        alter table org_api_keys
            add constraint org_api_keys_role_chk
            check (role in ('admin', 'member', 'viewer'));
    end if;
end $$;
