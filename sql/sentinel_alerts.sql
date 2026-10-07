-- Vibrion Sentinel: recurrence / cluster watcher alerts.
--
-- One row per rule firing, deduplicated by dedupe_key so a re-run of the
-- watcher never double-alerts on the same event. A cluster that gains a new
-- member produces a NEW dedupe_key (the sorted accession set changed), which
-- is the intended "cluster grew" signal -- not a duplicate.
--
-- Run once in the Supabase SQL editor (or psql against the project's
-- connection string), alongside sql/sentinel_runs.sql.

create table if not exists sentinel_alerts (
    id              bigint generated always as identity primary key,
    rule_id         text not null,              -- 'R1' | 'R2'
    region          text not null default 'haiti',
    dedupe_key      text not null unique,
    accessions      text[] not null,
    payload         jsonb not null,             -- full finding from watch_rules.py
    notified        boolean not null default false,
    fired_at        timestamptz not null default now(),
    created_at      timestamptz not null default now()
);

create index if not exists sentinel_alerts_rule_fired_idx
    on sentinel_alerts (rule_id, fired_at desc);
