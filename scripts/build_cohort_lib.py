"""Pure logic for scripts/build_cohort.py: sorting/limiting ENA rows into
a deterministic subset, and formatting the committed manifest text.
Kept separate from ena_search() (the actual HTTP call) and from main()'s
argparse/file-write wiring so this can be tested without network access.
"""


def get_field(header, row, name):
    """value at column `name` in `row`, given the TSV header list. ""
    if the column isn't present, or the row is shorter than expected
    (ENA has, in practice, returned short rows for some fields)."""
    i = header.index(name) if name in header else -1
    return row[i] if 0 <= i < len(row) else ""


def sort_and_limit(header, rows, limit=0):
    """Rows sorted by (collection_date, run_accession) for a
    deterministic, reproducible manifest -- ENA's own result order is not
    guaranteed stable across requests. limit=0 keeps every row."""
    rows = sorted(
        rows,
        key=lambda r: (get_field(header, r, "collection_date"),
                        get_field(header, r, "run_accession")),
    )
    if limit:
        rows = rows[:limit]
    return rows


def date_range(header, rows):
    """(min, max) collection_date across rows with a non-empty date, or
    None if none of them have one."""
    dates = [d for d in (get_field(header, r, "collection_date") for r in rows) if d]
    return (min(dates), max(dates)) if dates else None


def format_manifest(cohort_name, query, header, rows):
    """The exact committed-TSV text: a few '#'-prefixed metadata lines,
    then the header, then one row per line. Kept as a pure string
    builder so a test can assert on it directly rather than reading a
    file back."""
    lines = [
        f"# cohort: {cohort_name}",
        f"# ena_query: {query}",
        f"# runs: {len(rows)}",
        "\t".join(header),
    ]
    lines += ["\t".join(r) for r in rows]
    return "\n".join(lines) + "\n"
