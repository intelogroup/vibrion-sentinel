"""HL2601 fixture integrity: the deposited assembly stays what the
manifest says it is.

HL2601 (Front. Microbiol. 2026; GCF_059628145.1) is kept in this repo as
a VERIFICATION RECORD + FIXTURE -- never as a pipeline reference: no
production path accepts a deposited assembly (the assembly tier is
Nanopore-only, enforced in report_lib.validate_platform_config and
service validation), so this test deliberately does NOT feed the genome
through the pipeline expecting a classification. What it pins:

  * MANIFEST.tsv names the pinned accession (house manifest style).
  * The fixture fasta still measures exactly what the verification
    record states (total length, scaffold records, gap-split contigs).
  * The recorded contig count (51, per NCBI assembly stats) exceeds
    the repo assembly gate ASM_MAX_CONTIGS=50: this assembly would NOT
    clear the assembly-tier gate. Recorded, not retuned.
  * No assertion depends on PubMLST/tseemann database state (ST2040
    exists in PubMLST now; the oracle moves, the fixture must not).
    The recorded MLST row's tri-state reading belongs to the
    assembly-tier change and is tested there.
"""
import gzip
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "sentinel_lite"))
import report_lib  # noqa: E402

FIXTURES = Path(__file__).parent.parent / "fixtures" / "hl2601"


def fixture_stats():
    records, name, chunks = {}, None, []
    with gzip.open(FIXTURES / "genomic.fasta.gz", "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks).upper()
    total = sum(len(s) for s in records.values())
    contigs = 0
    for s in records.values():
        contigs += len([seg for seg in _split_gaps(s) if seg])
    lengths = sorted((len(s) for s in records.values()), reverse=True)
    acc, n50 = 0, 0
    for length in lengths:
        acc += length
        if acc >= total / 2:
            n50 = length
            break
    return {"records": len(records), "length_bp": total,
            "contigs": contigs, "n50_bp": n50}


def _split_gaps(seq):
    import re
    return re.split(r"N+", seq)


def manifest_row(filename):
    with open(FIXTURES / "MANIFEST.tsv") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0] == filename:
                return dict(zip(header, parts))
    return None


class TestHL2601Fixture(unittest.TestCase):
    def test_manifest_pins_accession_house_style(self):
        with open(FIXTURES / "MANIFEST.tsv") as f:
            header = f.readline().rstrip("\n")
        self.assertEqual(
            header, "file\taccession\tcoords\tlength_bp\tsource_paper\tverification")
        row = manifest_row("genomic.fasta.gz")
        self.assertIsNotNone(row)
        self.assertIn("GCF_059628145.1", row["accession"])
        self.assertIn("666", row["verification"])

    def test_fasta_measures_match_the_record(self):
        stats = fixture_stats()
        self.assertEqual(stats["length_bp"], 4297354)
        self.assertEqual(stats["records"], 47)   # scaffold records
        self.assertEqual(stats["contigs"], 51)   # split on N-gaps (NCBI)
        self.assertEqual(stats["n50_bp"], 324224)
        row = manifest_row("genomic.fasta.gz")
        self.assertEqual(int(row["length_bp"]), stats["length_bp"])

    def test_recorded_contig_count_would_fail_assembly_gate(self):
        # The honest gate fact, pinned so nobody "fixes" it later by
        # retuning ASM_MAX_CONTIGS: at NCBI's contig count of 51 this
        # assembly fails the gate on contigs alone (length and N50 are
        # inside the repo's constants: max length 4.4 Mb, N50 >= 200 kb).
        stats = fixture_stats()
        gate_stats = {"length_bp": stats["length_bp"],
                      "contigs": stats["contigs"],
                      "n50_bp": stats["n50_bp"]}
        self.assertLessEqual(gate_stats["length_bp"],
                             report_lib.ASM_MAX_LENGTH_BP)
        status, reasons = report_lib.assembly_qc_gate(gate_stats)
        self.assertEqual(status, "fail")
        self.assertTrue(any("contig" in r for r in reasons), reasons)
        self.assertEqual(report_lib.ASM_MAX_CONTIGS, 50)


if __name__ == "__main__":
    unittest.main()
