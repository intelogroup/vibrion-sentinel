#!/usr/bin/env python3
"""
Call variants (SNPs/indels) from aligned BAM file
Uses bcftools mpileup + call for variant detection
"""

import subprocess
import json
import gzip
import sys
from pathlib import Path

def parse_vcf_for_surveillance(vcf_path):
    """
    Parse VCF file and extract surveillance-relevant SNPs
    Focus on known resistance/virulence loci
    """
    
    # Surveillance loci of interest with their genomic coordinates
    # These correspond to genes we're watching for mutations
    SURVEILLANCE_GENES = {
        'gyrA': (2420000, 2423000, 'Fluoroquinolone resistance'),
        'parC': (965000, 968000, 'Fluoroquinolone resistance'),
        'rfb': (1330000, 1370000, 'O-antigen/vaccine escape'),
        'ctxAB': (437000, 445000, 'Cholera toxin'),
        'tcpA': (522000, 525000, 'Toxin co-regulated pilus'),
        'vps': (285000, 315000, 'Biofilm/rugose phenotype'),
        'hapR': (2650000, 2653000, 'Quorum sensing'),
        'rpsL': (2950000, 2953000, 'Streptomycin resistance'),
        'fusA': (2070000, 2074000, 'Fusidic acid resistance'),
    }
    
    snps = []
    functional_snps = []
    indels = []
    
    # Parse VCF
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    
    with open_func(vcf_path, 'rt') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
                
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            qual = float(fields[5]) if fields[5] != '.' else 0
            info = fields[7]
            
            # Skip low-quality variants (Snippy standard)
            if qual < 30:
                continue
            
            # Parse depth and allele frequency from INFO field
            depth = 0
            allele_freq = 0.0
            for item in info.split(';'):
                if item.startswith('DP='):
                    depth = int(item.split('=')[1])
                elif item.startswith('AF='):
                    allele_freq = float(item.split('=')[1].split(',')[0])
            
            # Apply Snippy-compatible filters for high-quality SNP calling
            # Min depth 10x, min AF 0.9 for clonal haploid organisms
            if depth < 10:
                continue
            if allele_freq < 0.9 and allele_freq > 0:  # Skip if AF=0 (not called)
                continue
            
            # Determine variant type
            is_snp = len(ref) == 1 and len(alt) == 1
            is_indel = len(ref) != len(alt)
            
            # Check if variant is in a surveillance locus
            gene_context = None
            
            # 🧬 REFERENCE AGNOSTIC COORDINATE MAPPING
            # O1 and O139 have different coordinates. 
            # If chrom starts with 'LC' (O139), use O139 mapping.
            is_o139 = chrom.startswith('LC594838')
            
            if is_o139:
                O139_MAPPING = {
                    'wbeT': (0, 38000, 'Serogroup O139 (wbf cluster)'), # O139 specific
                }
                for gene, (start, end, description) in O139_MAPPING.items():
                    if start <= pos <= end:
                        gene_context = {'gene': gene, 'description': description}
                        break
            else:
                for gene, (start, end, description) in SURVEILLANCE_GENES.items():
                    if start <= pos <= end:
                        gene_context = {
                            'gene': gene,
                            'description': description
                        }
                        break
            
            variant_info = {
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'depth': depth,
                'gene_context': gene_context
            }
            
            if is_snp:
                snps.append(variant_info)
                if gene_context:
                    functional_snps.append(variant_info)
            elif is_indel:
                indels.append(variant_info)
    
    # Phase 3: Clonal Filter & 37k SNP Alarm
    MAX_7PET_SNPS = 37000
    NORMAL_DRIFT_PER_YEAR = 4.4
    reference_year = 2010
    sample_year = 2022  # Default for current Haiti surveillance
    
    years_diff = sample_year - reference_year
    snp_count = len(snps)
    velocity = snp_count / years_diff if years_diff > 0 else 0
    
    alarm_37k = snp_count > MAX_7PET_SNPS
    # Anomaly if current drift is >5x the expected endemic velocity
    velocity_anomaly = velocity > (NORMAL_DRIFT_PER_YEAR * 5)
    
    status_msg = "NORMAL"
    if alarm_37k:
        status_msg = "LINEAGE_REPLACEMENT"
    elif velocity_anomaly:
        status_msg = "RAPID_EVOLUTION"

    return {
        'total_snps': snp_count,
        'functional_snps': len(functional_snps),
        'total_indels': len(indels),
        'snp_details': snps[:50],  # Top 50 SNPs
        'functional_snp_details': functional_snps,
        'indel_details': indels[:20],  # Top 20 indels
        'clonal_filter': {
            'status': status_msg,
            '37k_alarm': alarm_37k,
            'velocity': velocity,
            'velocity_anomaly': velocity_anomaly,
            'reference_year': reference_year,
            'sample_year': sample_year
        }
    }


def call_variants(bam_path, reference_path, output_vcf, threads=4):
    """
    Call variants using bcftools mpileup + call
    """
    
    # Index reference if not already done
    fai_path = Path(f"{reference_path}.fai")
    if not fai_path.exists():
        print(f"Indexing reference: {reference_path}")
        subprocess.run(['samtools', 'faidx', reference_path], check=True)
    
    # Run bcftools mpileup + call pipeline
    print(f"Calling variants from {bam_path}...")
    
    mpileup_cmd = [
        'bcftools', 'mpileup',
        '-Ou',  # Uncompressed BCF output
        '-f', reference_path,
        '--threads', str(threads),
        '-q', '30',  # Min mapping quality (Snippy standard)
        '-Q', '30',  # Min base quality (Snippy standard)
        '-d', '1000',  # Max depth per sample
        bam_path
    ]
    
    call_cmd = [
        'bcftools', 'call',
        '-mv',  # Multiallelic and variant sites only
        '-Oz',  # Compressed VCF output
        '-o', output_vcf,
        '--threads', str(threads),
        '--ploidy', '1'  # Haploid organism
    ]
    
    # Filter variants after calling to match Snippy quality standards
    # Require: min depth 10, min allele frequency 0.9 (haploid clonal)
    filter_expr = 'DP>=10 && AF>=0.9'
    
    # Run mpileup and pipe to call
    mpileup_proc = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    call_proc = subprocess.Popen(call_cmd, stdin=mpileup_proc.stdout, stderr=subprocess.PIPE)
    
    mpileup_proc.stdout.close()
    stdout, stderr = call_proc.communicate()
    
    if call_proc.returncode != 0:
        if stderr:
            print(f"BCFTOOLS CALL ERROR:\n{stderr.decode()}", file=sys.stderr)
        
        # Check mpileup error too
        _, mpileup_stderr = mpileup_proc.communicate()
        if mpileup_stderr:
             print(f"BCFTOOLS MPILEUP ERROR:\n{mpileup_stderr.decode()}", file=sys.stderr)

        raise RuntimeError(f"Variant calling failed with return code {call_proc.returncode}")
    
    # Filter VCF to match Snippy/published standards
    print(f"Filtering variants: {filter_expr}")
    filtered_vcf = output_vcf.replace('.vcf.gz', '.filtered.vcf.gz')
    filter_cmd = [
        'bcftools', 'view',
        '-i', filter_expr,
        '-Oz', '-o', filtered_vcf,
        output_vcf
    ]
    subprocess.run(filter_cmd, check=True)
    
    # Replace original with filtered version
    subprocess.run(['mv', filtered_vcf, output_vcf], check=True)
    
    # Index VCF
    print(f"Indexing VCF: {output_vcf}")
    subprocess.run(['bcftools', 'index', '-t', output_vcf], check=True)
    
    print(f"✅ Variants called and saved to {output_vcf}")


def calculate_snp_distance_to_haiti_refs(snp_count, reference_used):
    """
    Calculate SNP distance to Haiti 2010 and 2022 references
    """
    haiti_refs = {
        'data/references/2010EL-1786.fasta': 'Haiti 2010',
        'data/references/Haiti_2022_Resurgence.fasta': 'Haiti 2022'
    }
    
    distances = {}
    
    # The reference used has 0 distance to itself by definition
    ref_name = haiti_refs.get(reference_used, Path(reference_used).stem)
    
    if reference_used in haiti_refs:
        # We aligned to one of the Haiti references
        distances[haiti_refs[reference_used]] = snp_count
        
        # Add the other Haiti reference (we don't have its SNP distance without re-alignment)
        for ref_path, ref_label in haiti_refs.items():
            if ref_path != reference_used:
                distances[ref_label] = None  # Would need re-alignment
    else:
        # We aligned to a different reference (e.g., Malawi)
        # SNP distances to Haiti refs would need re-alignment
        distances['Haiti 2010'] = None
        distances['Haiti 2022'] = None
    
    return {
        'reference_used': ref_name,
        'snp_distance_to_reference': snp_count,
        'haiti_distances': distances,
        'note': 'Distances to other references require re-alignment'
    }


def main(snakemake):
    """
    Main function called by Snakemake
    """
    
    bam_path = snakemake.input.bam
    reference_path = snakemake.input.reference
    output_vcf = snakemake.output.vcf
    snp_report_path = snakemake.output.snp_report
    threads = snakemake.threads
    
    # Create output directory
    output_dir = Path(output_vcf).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Call variants
    call_variants(bam_path, reference_path, output_vcf, threads)
    
    # Parse VCF and extract surveillance-relevant variants
    print("Parsing VCF for surveillance loci...")
    snp_report = parse_vcf_for_surveillance(output_vcf)
    
    # Add SNP distance analysis
    snp_distance_info = calculate_snp_distance_to_haiti_refs(
        snp_report['total_snps'],
        reference_path
    )
    snp_report['snp_distance_analysis'] = snp_distance_info
    
    # Save SNP report
    with open(snp_report_path, 'w') as f:
        json.dump(snp_report, f, indent=2)
    
    print(f"✅ SNP report saved to {snp_report_path}")
    print(f"   Total SNPs: {snp_report['total_snps']}")
    print(f"   Functional SNPs (in surveillance loci): {snp_report['functional_snps']}")
    print(f"   Indels: {snp_report['total_indels']}")
    print(f"   SNP distance to {snp_distance_info['reference_used']}: {snp_distance_info['snp_distance_to_reference']}")


if __name__ == '__main__':
    main(snakemake)
