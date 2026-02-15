
import unittest
import sys
import os
from pathlib import Path

# Add scripts to path
sys.path.append(str(Path(__file__).parent.parent / "workflow" / "scripts"))

from predict_serotype import predict_serotype, analyze_toxin, check_environmental_markers

class TestPredictSerotype(unittest.TestCase):

    def test_predict_serotype_ogawa(self):
        # Empty variants = Wild Type = Ogawa
        wbeT_vars = []
        result = predict_serotype(wbeT_vars)
        self.assertEqual(result['prediction'], "Ogawa")
        self.assertIn("Vaccine Match", result['status'])

    def test_predict_serotype_inaba(self):
        # Any variant in wbeT triggers Inaba suspicion
        wbeT_vars = [{'pos': 100, 'ref': 'G', 'alt': 'A'}]
        result = predict_serotype(wbeT_vars)
        self.assertIn("Inaba", result['prediction'])
        self.assertIn("Vaccine Mismatch", result['status'])

    def test_analyze_toxin_classical(self):
        # No variants vs Reference (ctxB7) = Classical
        ctxB_vars = []
        result = analyze_toxin(ctxB_vars)
        self.assertIn("ctxB7", result['genotype'])
        self.assertIn("Hypervirulent", result['virulence'])

    def test_analyze_toxin_variant(self):
        ctxB_vars = [{'pos': 20, 'ref': 'C', 'alt': 'T'}]
        result = analyze_toxin(ctxB_vars)
        self.assertEqual(result['genotype'], "Variant / Atypical")

if __name__ == '__main__':
    unittest.main()
