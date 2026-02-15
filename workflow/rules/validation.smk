# Forensic Core Validation Rules
# Housekeeping gene checksum, CTX integration, and SXT assembly

# Rule: Housekeeping Gene Checksum Validation
rule validate_checksum:
    input:
        vcf=rules.snpeff_annotate.output.annotated_vcf,
        bed="data/references/housekeeping_genes.bed"
    output:
        report="{output_dir}/{{sample}}/09_consensus/qc_checksum.json".format(
            output_dir=config["output_dir"]
        )
    params:
        sample_id=lambda wildcards: wildcards.sample,
        permissive_flag="--permissive" if config.get("checksum_permissive", False) else ""
    conda:
        "../envs/analysis.yaml"
    log:
        "{output_dir}/{{sample}}/logs/validate_checksum.log".format(output_dir=config["output_dir"])
    shell:
        """
        python3 workflow/scripts/validate_checksum.py \
            --vcf {input.vcf} \
            --bed {input.bed} \
            --output {output.report} \
            --sample {params.sample_id} \
            {params.permissive_flag} \
            2>&1 | tee {log}
        """

# Rule: CTXφ Dual-Site Integration Detection
rule detect_ctx_integration:
    input:
        bam=rules.align_to_reference.output.bam,
        bai=rules.align_to_reference.output.bai,
        reference=select_reference
    output:
        report="{output_dir}/{{sample}}/09_consensus/ctx_integration.json".format(
            output_dir=config["output_dir"]
        )
    params:
        sample_id=lambda wildcards: wildcards.sample,
        min_depth=config.get("ctx_min_depth", 5)
    conda:
        "../envs/analysis.yaml"
    log:
        "{output_dir}/{{sample}}/logs/detect_ctx_integration.log".format(output_dir=config["output_dir"])
    shell:
        """
        python3 workflow/scripts/detect_ctx_integration.py \
            --bam {input.bam} \
            --reference {input.reference} \
            --output {output.report} \
            --sample {params.sample_id} \
            --min-depth {params.min_depth} \
            2>&1 | tee {log}
        """

# Rule: SXT Element Local Assembly
rule assemble_sxt:
    input:
        bam=rules.align_to_reference.output.bam,
        bai=rules.align_to_reference.output.bai,
        reference=select_reference
    output:
        report="{output_dir}/{{sample}}/09_consensus/sxt_assembly.json".format(
            output_dir=config["output_dir"]
        ),
        contigs="{output_dir}/{{sample}}/09_consensus/sxt_contigs.fasta".format(
            output_dir=config["output_dir"]
        )
    params:
        sample_id=lambda wildcards: wildcards.sample,
        mode=config.get("pipeline_mode", "LABORATORY_FULL"),
        outdir=lambda wildcards: f"{config['output_dir']}/{wildcards.sample}/09_consensus/sxt_work",
        threads=config.get("threads", 4)
    conda:
        "../envs/assembly.yaml"
    log:
        "{output_dir}/{{sample}}/logs/assemble_sxt.log".format(output_dir=config["output_dir"])
    shell:
        """
        python3 workflow/scripts/assemble_sxt.py \
            --bam {input.bam} \
            --reference {input.reference} \
            --output {output.report} \
            --contigs {output.contigs} \
            --outdir {params.outdir} \
            --sample {params.sample_id} \
            --mode {params.mode} \
            --threads {params.threads} \
            2>&1 | tee {log}
        """
