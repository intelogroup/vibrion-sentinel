
import unittest
import sys
from pathlib import Path

# Add scripts to path
sys.path.append(str(Path(__file__).parent.parent / "workflow" / "scripts"))

from predict_serotype import predict_serotype, analyze_toxin, assess_hapR_status

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

    # === hapR Regulator Status Tests ===

    def test_hapR_functional(self):
        """No variants = FUNCTIONAL hapR (normal QS regulation)"""
        result = assess_hapR_status([], ctxAB_detected=False)
        self.assertEqual(result['regulator_status'], "FUNCTIONAL")
        self.assertEqual(result['threat_multiplier'], 1.0)

    def test_hapR_variant_snp(self):
        """SNP (same length ref/alt) = VARIANT (partial disruption)"""
        variants = [{'pos': 100, 'ref': 'G', 'alt': 'A', 'info': 'missense_variant'}]
        result = assess_hapR_status(variants, ctxAB_detected=False)
        self.assertEqual(result['regulator_status'], "VARIANT")
        self.assertEqual(result['threat_multiplier'], 1.1)

    def test_hapR_lof_without_ctxAB(self):
        """Frameshift (indel) without ctxAB = LOF_MUTANT (biofilm derepressed)"""
        variants = [{'pos': 200, 'ref': 'GC', 'alt': 'G', 'info': 'frameshift_variant'}]
        result = assess_hapR_status(variants, ctxAB_detected=False)
        self.assertEqual(result['regulator_status'], "LOF_MUTANT")
        self.assertEqual(result['threat_multiplier'], 1.3)

    def test_hapR_lof_with_ctxAB(self):
        """Frameshift + ctxAB = DEREPRESSED_VIRULENCE (super-producer)"""
        variants = [{'pos': 200, 'ref': 'GC', 'alt': 'G', 'info': 'frameshift_variant'}]
        result = assess_hapR_status(variants, ctxAB_detected=True)
        self.assertEqual(result['regulator_status'], "DEREPRESSED_VIRULENCE")
        self.assertEqual(result['threat_multiplier'], 1.5)
        self.assertIn("Super-producer", result['implication'])

    def test_hapR_stop_gained_with_ctxAB(self):
        """stop_gained in INFO + ctxAB = DEREPRESSED_VIRULENCE"""
        variants = [{'pos': 150, 'ref': 'C', 'alt': 'T', 'info': 'stop_gained'}]
        result = assess_hapR_status(variants, ctxAB_detected=True)
        self.assertEqual(result['regulator_status'], "DEREPRESSED_VIRULENCE")

if __name__ == '__main__':
    unittest.main()
