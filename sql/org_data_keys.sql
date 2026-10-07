-- Vibrion Sentinel: per-org data-encryption keys (Phase 2, envelope).
--
-- Key hierarchy: SENTINEL_KEK (env, AES-256) wraps each org's DEK (AES-256,
-- random per org). Upload bytes are encrypted with the DEK under AES-256-GCM
-- with a fresh 96-bit nonce per file; {dek_version, nonce} travel in a
-- per-object sidecar ({key}.meta.json), never in the database.
--
-- org_data_keys: the CURRENT DEK pointer per org.
-- org_data_key_versions: every DEK version ever issued. Old versions are
--   RETAINED (never deleted): history decrypts with its recorded version.
--   Rotation creates a new version; it does not re-encrypt history.
--
-- dek_wrapped is the DEK wrapped by the KEK: nonce(12) || ciphertext+tag,
-- stored base64-encoded. (Column type is bytea per the design; the REST
-- client sends base64, the postgrest-py convention.)
--
-- Run once in the Supabase SQL editor (or psql), after the Phase 1 tables.

create table if not exists org_data_key_versions (
    org_id      text not null,
    dek_version int not null,
    dek_wrapped bytea not null,   -- KEK-wrapped DEK (base64 over the REST API)
    kek_id      text not null,    -- which KEK wrapped it; 'env-v1' for now
    created_at  timestamptz not null default now(),
    primary key (org_id, dek_version)
);

create table if not exists org_data_keys (
    org_id      text primary key,
    dek_version int not null,
    dek_wrapped bytea not null,
    kek_id      text not null,
    created_at  timestamptz not null default now(),
    rotated_at  timestamptz,      -- set when this row stops being current
    foreign key (org_id, dek_version)
        references org_data_key_versions (org_id, dek_version)
);

create index if not exists org_data_key_versions_org_idx
    on org_data_key_versions (org_id);

-- Service-role writes only (service_role bypasses RLS anyway); read policy
-- mirrors the other service tables.
alter table org_data_key_versions enable row level security;
alter table org_data_keys enable row level security;

drop policy if exists "org_data_key_versions read access" on org_data_key_versions;
create policy "org_data_key_versions read access"
    on org_data_key_versions for select using (true);

drop policy if exists "org_data_keys read access" on org_data_keys;
create policy "org_data_keys read access"
    on org_data_keys for select using (true);
