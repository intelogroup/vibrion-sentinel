-- Vibrion Sentinel: user sessions (Phase 3 portal).
--
-- Opaque tokens: sha256 hex at rest, plaintext only in the HttpOnly cookie
-- and the one-time login/signup JSON response. 30-day expiry; revoked on
-- logout and on password change (all sessions die then).
--
-- Run this once in the Supabase SQL editor (or `psql`).

create table if not exists user_sessions (
    id          uuid primary key default gen_random_uuid(),
    user_id     uuid not null references users (id) on delete cascade,
    token_hash  text not null unique,        -- sha256 hex of the opaque token
    created_at  timestamptz not null default now(),
    expires_at  timestamptz not null,
    revoked_at  timestamptz                  -- null = active
);

create index if not exists user_sessions_token_hash_idx on user_sessions (token_hash);
create index if not exists user_sessions_user_id_idx on user_sessions (user_id);

alter table user_sessions enable row level security;

drop policy if exists "user_sessions read access" on user_sessions;
create policy "user_sessions read access"
    on user_sessions for select
    using (true);
