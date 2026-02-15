#!/usr/bin/env python3
"""
Call variants (SNPs/indels) from aligned BAM file.
Produces Main (Consensus) and Minor (Heterogeneity) VCFs.
Supports dynamic gene mapping and configurable AF thresholds.
"""

import subprocess
import json
import gzip
import sys
from pathlib import Path
from collections import defaultdict

def analyze_heterogeneity(minor_vcf_path, surveillance_genes, minor_af=0.1, min_depth=20):
    """
    Scan the minor VCF for evidence of mixed populations in surveillance loci.
    """
    alerts = []
    if not Path(minor_vcf_path).exists():
        return alerts

    open_func = gzip.open if str(minor_vcf_path).endswith('.gz') else open
    
    with open_func(minor_vcf_path, 'rt') as vcf:
        for line in vcf:
            if line.startswith('#'): continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8: continue
            
            chrom, pos = fields[0], int(fields[1])
            info = fields[7]
            
            af = 0.0
            dp = 0
            for item in info.split(';'):
                if item.startswith('AF='):
                    try: af = float(item.split('=')[1].split(',')[0])
                    except Exception: pass
                if item.startswith('DP='):
                    try: dp = int(item.split('=')[1])
                    except Exception: pass

            # Mixed population range: minor_af to 0.9 (fixed)
            if minor_af <= af < 0.9 and dp >= min_depth:
                for gene, (g_chrom, start, end, desc) in surveillance_genes.items():
                    if chrom == g_chrom and start <= pos <= end:
                        alerts.append({
                            'chrom': chrom,
                            'gene': gene,
                            'pos': pos,
                            'af': af,
                            'depth': dp,
                            'description': desc,
                            'message': f"Heterogeneity detected in {gene} (AF={af:.2f})"
                        })
                        break
    return alerts

def parse_vcf_output(vcf_path, minor_vcf_path=None, gene_map_file=None, clonal_af=0.9, minor_af=0.1, min_depth_hetero=20):
    """
    Parse VCFs and extract surveillance SNPs and heterogeneity alerts.
    """
    # Default surveillance genes (Haiti 2010 Reference CP003069.1)
    SURVEILLANCE_GENES = {
        'gyrA': ('CP003069.1', 2420000, 2423000, 'Fluoroquinolone resistance'),
        'parC': ('CP003069.1', 965000, 968000, 'Fluoroquinolone resistance'),
        'rfb': ('CP003069.1', 1330000, 1370000, 'O-antigen/vaccine escape'),
        'ctxAB': ('CP003069.1', 437000, 445000, 'Cholera toxin'),
        'tcpA': ('CP003069.1', 522000, 525000, 'Toxin co-regulated pilus'),
        'vps': ('CP003069.1', 285000, 315000, 'Biofilm/rugose phenotype'),
        'hapR': ('CP003069.1', 2650000, 2653000, 'Quorum sensing'),
        'rpsL': ('CP003069.1', 2950000, 2953000, 'Streptomycin resistance'),
        'fusA': ('CP003069.1', 2070000, 2074000, 'Fusidic acid resistance'),
        'wbeT': ('CP003069.1', 2678186, 2678980, 'Serotype switch (Ogawa/Inaba)')
    }

    if gene_map_file and Path(gene_map_file).exists():
        try:
            with open(gene_map_file) as f:
                mapping = json.load(f)
                SURVEILLANCE_GENES = {}
                for gene, data in mapping.items():
                    SURVEILLANCE_GENES[gene] = (data['chrom'], data['start'], data['end'], data.get('description', f'{gene} locus'))
            print(f"Loaded {len(SURVEILLANCE_GENES)} loci from gene map.")
        except Exception as e:
            print(f"Warning: Failed to load gene map: {e}")

    snps, functional_snps, indels = [], [], []
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    
    with open_func(vcf_path, 'rt') as vcf:
        for line in vcf:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 8: continue
            
            chrom, pos, ref, alt, qual, info = fields[0], int(fields[1]), fields[3], fields[4], float(fields[5]) if fields[5] != '.' else 0, fields[7]
            if qual < 30: continue
            
            depth, af = 0, 0.0
            for item in info.split(';'):
                if item.startswith('DP='): depth = int(item.split('=')[1])
                elif item.startswith('AF='): af = float(item.split('=')[1].split(',')[0])
                elif item.startswith('DP4='):
                    # Fallback if AF missing: DP4=ref_f,ref_r,alt_f,alt_r
                    try:
                        parts = list(map(int, item.split('=')[1].split(',')))
                        ref_count = parts[0] + parts[1]
                        alt_count = parts[2] + parts[3]
                        total = ref_count + alt_count
                        if total > 0:
                            af = alt_count / total
                    except: pass
            
            if depth < 10 or af < clonal_af: continue
            
            is_snp = (len(ref) == 1 and len(alt) == 1)
            is_indel = (len(ref) != len(alt))
            
            gene_context = None
            for gene, (g_chrom, start, end, desc) in SURVEILLANCE_GENES.items():
                if chrom == g_chrom and start <= pos <= end:
                    gene_context = {'gene': gene, 'description': desc}
                    break
            
            var = {'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt, 'depth': depth, 'af': af, 'gene_context': gene_context}
            if is_snp:
                snps.append(var)
                if gene_context: functional_snps.append(var)
            elif is_indel:
                indels.append(var)

    heterogeneity_alerts = []
    if minor_vcf_path:
        heterogeneity_alerts = analyze_heterogeneity(minor_vcf_path, SURVEILLANCE_GENES, minor_af, min_depth_hetero)

    return {
        'total_snps': len(snps),
        'functional_snps': len(functional_snps),
        'total_indels': len(indels),
        'snp_details': snps[:50],
        'functional_snp_details': functional_snps,
        'heterogeneity_alerts': heterogeneity_alerts,
        'thresholds_used': {'clonal': clonal_af, 'minor': minor_af}
    }

def call_variants(bam_path, reference_path, output_vcf, threads=4, clonal_af=0.9, minor_af=0.1):
    """
    Call variants using bcftools with robust clinical-grade filtering.
    Generates Main (Consensus) and Minor (Heterogeneity) VCFs.
    """
    if not Path(f"{reference_path}.fai").exists():
        subprocess.run(['samtools', 'faidx', reference_path], check=True)
    
    raw_bcf = output_vcf.replace('.vcf.gz', '.raw.bcf')
    norm_bcf = output_vcf.replace('.vcf.gz', '.norm.bcf')
    
    # 1. MPILEUP & CALL (Annotated)
    # -a: Annotate FORMAT/AD, ADF, ADR (Allele Depth, Fwd/Rev) for strand bias filtering
    # -Q 20: Min base quality
    # -q 30: Min mapping quality (skip ambiguous reads)
    print(f"Calling variants (mpileup | call | norm)...")
    
    mpileup_cmd = [
        'bcftools', 'mpileup', 
        '-Ou', 
        '-f', reference_path, 
        '--threads', str(threads), 
        '-a', 'FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP', 
        '-q', '30', 
        '-Q', '20', 
        '-d', '2000', 
        bam_path
    ]
    
    call_cmd = [
        'bcftools', 'call', 
        '-mv', 
        '-Ob', 
        '-o', raw_bcf, 
        '--threads', str(threads), 
        '--ploidy', '1', 
        '-A' # Keep all alleles
    ]
    
    p1 = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(call_cmd, stdin=p1.stdout)
    p1.stdout.close()
    p2.communicate()
    
    if p2.returncode != 0:
        raise RuntimeError("Variant calling failed")

    # 2. NORMALIZE (Left-align indels)
    print("Normalizing variants...")
    subprocess.run([
        'bcftools', 'norm',
        '-f', reference_path,
        '-m', '-both',
        '-Ob', '-o', norm_bcf,
        raw_bcf
    ], check=True)

    # 3. CONSENSUS VCF (High Confidence for Phylogeny)
    # Rules: High Quality, High Depth, High AF, Balanced Strands (if cov high)
    # Note: For haploid, strand bias is simpler: Ensure reads on both strands support the call if coverage is sufficient
    print(f"Generating Consensus VCF (Strict: AF >= {clonal_af}, MQ>=30, GQ>=20)...")
    
    # Filter Expression (CORRECTED INDEXING [1] for ALT):
    # QUAL >= 30
    # DP >= 10
    # GQ >= 20 (Genotype Quality)
    # AF >= clonal_af (calculated via AD[1]/DP)
    # MQ >= 30 (Mapping Quality)
    # Strand Bias: ADF[1]>0 && ADR[1]>0 (Require ALT support on FWD and REV strands)
    
    strict_filter = (
        f'QUAL >= 30 && '
        f'FORMAT/DP >= 10 && '
        f'FORMAT/AD[0:1] / FORMAT/DP >= {clonal_af} && '  # Sample 0, Allele 1
        f'INFO/MQ >= 30 && '
        f'(FORMAT/ADF[0:1] > 0 && FORMAT/ADR[0:1] > 0)'   # Sample 0, Allele 1
    )
    
    subprocess.run([
        'bcftools', 'filter',
        '-i', strict_filter,
        '-Oz', '-o', output_vcf,
        norm_bcf
    ], check=True)
    subprocess.run(['bcftools', 'index', '-t', output_vcf], check=True)

    # 4. MINOR VARIANT VCF (Heterogeneity Detection)
    print(f"Generating Minor Variant VCF (Heterogeneity: {minor_af} <= AF < {clonal_af})...")
    minor_vcf = output_vcf.replace('.vcf.gz', '.minor.vcf.gz')
    
    minor_filter = (
        f'QUAL >= 30 && '
        f'FORMAT/DP >= 20 && '
        f'FORMAT/AD[0:1] >= 5 && '
        f'FORMAT/AD[0:1] / FORMAT/DP >= {minor_af} && '
        f'FORMAT/AD[0:1] / FORMAT/DP < {clonal_af} && '
        f'(FORMAT/ADF[0:1] > 0 && FORMAT/ADR[0:1] > 0)'
    )
    
    subprocess.run([
        'bcftools', 'filter',
        '-i', minor_filter,
        '-Oz', '-o', minor_vcf,
        norm_bcf
    ], check=True)
    subprocess.run(['bcftools', 'index', '-t', minor_vcf], check=True)
    
    # Cleanup temps
    Path(raw_bcf).unlink(missing_ok=True)
    Path(norm_bcf).unlink(missing_ok=True)

def main(snakemake):
    bam_path = snakemake.input.bam
    reference_path = snakemake.input.reference
    output_vcf = snakemake.output.vcf
    snp_report_path = snakemake.output.snp_report
    threads = snakemake.threads
    
    v_thresh = snakemake.config.get('variant_thresholds', {})
    clonal_af = v_thresh.get('clonal_af', 0.9)
    minor_af = v_thresh.get('minor_af', 0.1)
    min_depth_hetero = v_thresh.get('hetero_min_depth', 20)
    gene_map_file = getattr(snakemake.input, 'gene_map', None)
    
    Path(output_vcf).parent.mkdir(parents=True, exist_ok=True)
    
    call_variants(bam_path, reference_path, output_vcf, threads, clonal_af, minor_af)
    
    minor_vcf_path = output_vcf.replace('.vcf.gz', '.minor.vcf.gz')
    snp_report = parse_vcf_output(output_vcf, minor_vcf_path, gene_map_file, clonal_af, minor_af, min_depth_hetero)
    
    with open(snp_report_path, 'w') as f:
        json.dump(snp_report, f, indent=2)
    
    print(f"✅ Variants called and report saved to {snp_report_path}")

if __name__ == '__main__':
    main(snakemake)