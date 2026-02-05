#!/usr/bin/env python3
"""
Calculate SNP distances to multiple Haiti references
Performs quick re-alignment and variant calling against Haiti 2010 and 2022
"""

import subprocess
import tempfile
import json
from pathlib import Path

def count_snps_from_vcf(vcf_path):
    """Count total SNPs in a VCF file"""
    result = subprocess.run(
        f'bcftools view -v snps {vcf_path} | bcftools query -f "%CHROM\\t%POS\\n" | wc -l',
        shell=True,
        capture_output=True,
        text=True
    )
    count_str = result.stdout.strip()
    return int(count_str) if count_str else 0


def quick_align_and_call(reads_fastq, reference_fasta, threads=4):
    """
    Quick alignment and variant calling
    Returns number of SNPs
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_bam = Path(tmpdir) / "tmp.bam"
        tmp_vcf = Path(tmpdir) / "tmp.vcf.gz"
        
        # Check if reference is indexed
        if not Path(f"{reference_fasta}.bwt").exists():
            print(f"   Indexing {reference_fasta}...")
            subprocess.run(['bwa', 'index', reference_fasta], 
                          check=True, capture_output=True)
        
        # Align
        print(f"   Aligning to {Path(reference_fasta).stem}...")
        with open(tmp_bam, 'wb') as bam_out:
            bwa = subprocess.Popen(
                ['bwa', 'mem', '-t', str(threads), reference_fasta, reads_fastq],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL
            )
            samtools = subprocess.Popen(
                ['samtools', 'view', '-b', '-'],
                stdin=bwa.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL
            )
            sort = subprocess.Popen(
                ['samtools', 'sort', '-@', str(threads), '-o', str(tmp_bam), '-'],
                stdin=samtools.stdout,
                stderr=subprocess.DEVNULL
            )
            sort.wait()
        
        # Index BAM
        subprocess.run(['samtools', 'index', str(tmp_bam)], 
                      check=True, capture_output=True)
        
        # Call variants with Snippy-compatible parameters
        print("   Calling variants (Q30, depth≥10, AF≥0.9)...")
        tmp_vcf_raw = Path(tmpdir) / "tmp_raw.vcf.gz"
        
        with open(tmp_vcf_raw, 'wb') as vcf_out:
            mpileup = subprocess.Popen(
                ['bcftools', 'mpileup', '-Ou', '-f', reference_fasta,
                 '-q', '30',  # Min mapping quality (Snippy standard)
                 '-Q', '30',  # Min base quality (Snippy standard)
                 '-d', '1000',  # Max depth per sample
                 '--ploidy', '1',  # Haploid organism
                 str(tmp_bam)],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL
            )
            call = subprocess.Popen(
                ['bcftools', 'call', '-mv', '--ploidy', '1', '-Oz', '-o', str(tmp_vcf_raw)],
                stdin=mpileup.stdout,
                stderr=subprocess.DEVNULL
            )
            call.wait()
        
        # Filter variants: DP>=10 && AF>=0.9
        filter_expr = 'DP>=10 && AF>=0.9'
        subprocess.run(
            ['bcftools', 'view', '-i', filter_expr, '-Oz', '-o', str(tmp_vcf), str(tmp_vcf_raw)],
            check=True,
            capture_output=True
        )
        
        # Index VCF
        subprocess.run(['bcftools', 'index', '-t', str(tmp_vcf)], 
                      check=True, capture_output=True)
        
        # Count SNPs
        snp_count = count_snps_from_vcf(str(tmp_vcf))
        
        return snp_count


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Calculate SNP distances to Haiti references')
    parser.add_argument('--reads', required=True, help='Input FASTQ file (Vibrio reads)')
    parser.add_argument('--output', required=True, help='Output JSON file')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('--sample', required=True, help='Sample name')
    
    args = parser.parse_args()
    
    haiti_refs = {
        'Haiti_2010': 'data/references/2010EL-1786.fasta',
        'Haiti_2022': 'data/references/Haiti_2022_Resurgence.fasta'
    }
    
    print(f"🧬 Calculating SNP distances for {args.sample}...")
    
    distances = {}
    
    for ref_name, ref_path in haiti_refs.items():
        if not Path(ref_path).exists():
            print(f"   ⚠️  {ref_name} reference not found: {ref_path}")
            distances[ref_name] = None
            continue
        
        print(f"\n📊 Testing {ref_name}...")
        snp_count = quick_align_and_call(args.reads, ref_path, args.threads)
        distances[ref_name] = snp_count
        print(f"   ✅ {snp_count} SNPs")
    
    # Save results
    output = {
        'sample': args.sample,
        'snp_distances': distances,
        'note': 'SNP counts relative to each reference genome',
        'filters': {
            'quality': 'Q30',
            'min_depth': 10,
            'min_allele_freq': 0.9,
            'max_depth': 1000,
            'method': 'Snippy-compatible'
        }
    }
    
    with open(args.output, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n✅ SNP distance report saved to {args.output}")
    print("\n📊 Summary:")
    for ref, dist in distances.items():
        if dist is not None:
            print(f"   {ref}: {dist} SNPs")
        else:
            print(f"   {ref}: N/A")


if __name__ == '__main__':
    main()
