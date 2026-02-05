# Consensus Polishing Rules (Platform-Aware)
# Insert these after the generate_consensus rule in Snakefile

# Rule 7b: Platform Detection (Auto-detect Nanopore vs Illumina)
rule detect_platform:
    input:
        fastq=rules.hostile_clean.output.cleaned
    output:
        report="{output_dir}/{{sample}}/09_consensus/platform_detection.json".format(
            output_dir=config["output_dir"]
        )
    conda:
        "../envs/analysis.yaml"
    log:
        "{output_dir}/{{sample}}/logs/platform_detection.log".format(output_dir=config["output_dir"])
    script:
        "../scripts/detect_platform.py"

# Rule 7c: Polish Consensus (Platform-Aware)
# Uses Medaka for Nanopore, Pilon for Illumina
rule polish_consensus:
    input:
        draft=rules.generate_consensus.output.fasta,
        platform_report=rules.detect_platform.output.report,
        reads=rules.hostile_clean.output.cleaned,
        bam=rules.align_to_reference.output.bam,
        coverage_report=rules.verify_coverage_integrity.output.report
    output:
        polished="{output_dir}/{{sample}}/09_consensus/{{sample}}_polished.fasta".format(
            output_dir=config["output_dir"]
        ),
        manifest="{output_dir}/{{sample}}/09_consensus/polishing_manifest.json".format(
            output_dir=config["output_dir"]
        )
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.polished),
        mode=config.get("pipeline_mode", "LABORATORY_FULL")
    conda:
        "../envs/analysis.yaml"
    threads: 8
    log:
        "{output_dir}/{{sample}}/logs/polish_consensus.log".format(output_dir=config["output_dir"])

    shell:
        """
        # Coverage Logic Gate (Thermal + Forensic Safety)
        COVERAGE=$(python3 -c "import json; print(json.load(open('{input.coverage_report}'))['metrics']['global_depth'])")
        MIN_COV=10.0
        
        echo "Coverage Check: ${{COVERAGE}}x (threshold: ${{MIN_COV}}x)" >> {log}
        
        # Check if coverage is sufficient for polishing
        if python3 -c "import sys; sys.exit(0 if $COVERAGE < $MIN_COV else 1)"; then
            echo "⚠️  Coverage too low for polishing (<${{MIN_COV}}x)" >> {log}
            echo "ℹ️  Thermal Safety: Skipping expensive polishing on low-coverage data" >> {log}
            echo "ℹ️  Forensic Safety: Preventing artifactual 'polished' genome from insufficient data" >> {log}
            echo "STATUS: LOW_COV_SKIP" >> {log}
            
            # Symlink draft to polished (pipeline continues)
            cp {input.draft} {output.polished}
            echo "✅ Draft consensus copied as final (no polishing applied)" >> {log}
        else
            echo "✅ Coverage sufficient (${{COVERAGE}}x >= ${{MIN_COV}}x). Proceeding with polishing." >> {log}
            
            # Read platform from detection report
            PLATFORM=$(python3 -c "import json; print(json.load(open('{input.platform_report}'))['platform'])")
            POLISHER=$(python3 -c "import json; print(json.load(open('{input.platform_report}'))['polisher_config']['tool'])")
            
            echo "Platform: $PLATFORM" >> {log}
            echo "Polisher: $POLISHER" >> {log}
            
            if [ "$POLISHER" = "medaka" ]; then
                # Nanopore polishing
                python3 workflow/scripts/polish_medaka.py \
                    --draft {input.draft} \
                    --reads {input.reads} \
                    --output {output.polished} \
                    --outdir {params.outdir}/medaka_work \
                    --rounds 2 \
                    2>&1 | tee -a {log}
            elif [ "$POLISHER" = "pilon" ]; then
                # Illumina polishing (Hybrid v2.1: Polypolish -> Pilon)
                python3 workflow/scripts/polish_hybrid_illumina.py \
                    --draft {input.draft} \
                    --reads {input.reads} \
                    --output {output.polished} \
                    --outdir {params.outdir}/hybrid_work \
                    --threads {threads} \
                    2>&1 | tee -a {log}
            else
                echo "Unknown polisher: $POLISHER" >> {log}
                cp {input.draft} {output.polished}
            fi
        fi

        # Export manifest variables
        export MANIFEST_INPUT_DRAFT="{input.draft}"
        export MANIFEST_INPUT_READS="{input.reads}"
        export MANIFEST_INPUT_BAM="{input.bam}"
        export MANIFEST_INPUT_PLATFORM="{input.platform_report}"
        export MANIFEST_INPUT_COVERAGE="{input.coverage_report}"
        export MANIFEST_OUTPUT_POLISHED="{output.polished}"
        export MANIFEST_OUTPUT_PATH="{output.manifest}"
        export MANIFEST_LOG="{log}"

        python3 - <<'PY'
import json
import hashlib
import subprocess
import time
import os

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def tool_version(cmd):
    try:
        output = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
        return output.splitlines()[0].strip()
    except Exception:
        return None

# Load coverage and platform reports
with open(os.environ['MANIFEST_INPUT_COVERAGE']) as f:
    coverage = json.load(f)

with open(os.environ['MANIFEST_INPUT_PLATFORM']) as f:
    platform_info = json.load(f)

manifest = {{
    "step": "polishing",
    "timestamp": time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    "inputs": {{
        "draft": os.environ['MANIFEST_INPUT_DRAFT'],
        "reads": os.environ['MANIFEST_INPUT_READS'],
        "bam": os.environ['MANIFEST_INPUT_BAM'],
        "platform_report": os.environ['MANIFEST_INPUT_PLATFORM'],
        "coverage_report": os.environ['MANIFEST_INPUT_COVERAGE']
    }},
    "outputs": {{
        "polished": os.environ['MANIFEST_OUTPUT_POLISHED']
    }},
    "hashes": {{
        "draft_sha256": sha256(os.environ['MANIFEST_INPUT_DRAFT']),
        "reads_sha256": sha256(os.environ['MANIFEST_INPUT_READS']),
        "bam_sha256": sha256(os.environ['MANIFEST_INPUT_BAM']),
        "platform_report_sha256": sha256(os.environ['MANIFEST_INPUT_PLATFORM']),
        "coverage_report_sha256": sha256(os.environ['MANIFEST_INPUT_COVERAGE']),
        "polished_sha256": sha256(os.environ['MANIFEST_OUTPUT_POLISHED'])
    }},
    "qc_metrics": {{
        "platform": platform_info.get("platform"),
        "polisher": platform_info.get("polisher_config", {{}}).get("tool"),
        "global_depth": coverage.get("metrics", {{}}).get("global_depth"),
        "wbeT_coverage": coverage.get("metrics", {{}}).get("wbeT_coverage"),
        "wbeT_mean_depth": coverage.get("metrics", {{}}).get("wbeT_mean_depth"),
        "polishing_applied": coverage.get("metrics", {{}}).get("global_depth", 0) >= 10.0
    }},
    "tools": {{
        "medaka": tool_version(["medaka", "--version"]),
        "pilon": tool_version(["pilon", "--version"]),
        "polypolish": tool_version(["polypolish", "--version"])
    }},
    "log": os.environ['MANIFEST_LOG']
}}

with open(os.environ['MANIFEST_OUTPUT_PATH'], "w") as f:
    json.dump(manifest, f, indent=2)
PY
        """
