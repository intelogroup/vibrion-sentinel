
import unittest
import sys
import tempfile
import gzip
from pathlib import Path

# Add scripts to path
sys.path.append(str(Path(__file__).parent.parent / "workflow" / "scripts"))

from call_variants import analyze_heterogeneity

class TestCallVariants(unittest.TestCase):

    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.vcf_path = Path(self.test_dir.name) / "test.minor.vcf.gz"
        
        # Create a mock VCF with a mixed variant
        # Gene: wbeT is around 2678186-2678980
        # Variant at 2678200 with AF=0.4, DP=100
        content = [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "CP003069.1\t2678200\t.\tG\tA\t999\t.\tDP=100;AF=0.4"
        ]
        
        with gzip.open(self.vcf_path, 'wt') as f:
            f.write("\n".join(content) + "\n")
            
        self.surveillance_genes = {
            'wbeT': ('CP003069.1', 2678186, 2678980, 'Serotype switch')
        }

    def tearDown(self):
        self.test_dir.cleanup()

    def test_heterogeneity_detection(self):
        alerts = analyze_heterogeneity(self.vcf_path, self.surveillance_genes, minor_af=0.1, min_depth=20)
        self.assertEqual(len(alerts), 1)
        self.assertEqual(alerts[0]['gene'], 'wbeT')
        self.assertAlmostEqual(alerts[0]['af'], 0.4)
        self.assertEqual(alerts[0]['depth'], 100)

    def test_heterogeneity_low_depth(self):
        # Update VCF with low depth
        content = [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "CP003069.1\t2678200\t.\tG\tA\t999\t.\tDP=10;AF=0.4"
        ]
        with gzip.open(self.vcf_path, 'wt') as f:
            f.write("\n".join(content) + "\n")
            
        alerts = analyze_heterogeneity(self.vcf_path, self.surveillance_genes, minor_af=0.1, min_depth=20)
        self.assertEqual(len(alerts), 0) # Should filter out

if __name__ == '__main__':
    unittest.main()
