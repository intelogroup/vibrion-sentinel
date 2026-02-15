
import unittest
import sys
from pathlib import Path
from collections import Counter

# Add scripts to path
sys.path.append(str(Path(__file__).parent.parent / "workflow" / "scripts"))

from resolve_haplotypes import get_consensus_from_reads

class MockRead:
    def __init__(self, sequence, start_pos):
        self.query_sequence = sequence
        self.reference_start = start_pos
        
    def get_aligned_pairs(self, matches_only=True):
        # Simple 1-to-1 mapping for test
        return [(i, self.reference_start + i) for i in range(len(self.query_sequence))]

class TestResolveHaplotypes(unittest.TestCase):

    def test_get_consensus_perfect(self):
        reads = [
            MockRead("ATCG", 100),
            MockRead("ATCG", 100),
            MockRead("ATCG", 100)
        ]
        consensus = get_consensus_from_reads(reads, 100, 104)
        self.assertEqual(consensus, "ATCG")

    def test_get_consensus_mixed(self):
        reads = [
            MockRead("ATCG", 100), # Major
            MockRead("ATCG", 100), # Major
            MockRead("ATCA", 100)  # Minor
        ]
        consensus = get_consensus_from_reads(reads, 100, 104)
        self.assertEqual(consensus, "ATCG") # Majority wins

    def test_get_consensus_gaps(self):
        reads = [
            MockRead("A", 100),
            MockRead("G", 103)
        ]
        # 100: A, 101: N, 102: N, 103: G
        consensus = get_consensus_from_reads(reads, 100, 104)
        self.assertEqual(consensus, "ANNG")

if __name__ == '__main__':
    unittest.main()
