#!/usr/bin/env python3
"""
Automated Polishing Integration Script
Safely integrates polishing rules into the Vibrion Snakefile
"""

import re
import sys
from pathlib import Path


def backup_file(filepath: Path) -> Path:
    """Create a backup of the file."""
    backup = filepath.with_suffix(filepath.suffix + '.backup')
    backup.write_text(filepath.read_text())
    print(f"✅ Backup created: {backup}")
    return backup


def insert_polishing_rules(snakefile_path: Path, rules_path: Path) -> bool:
    """Insert polishing rules after generate_consensus."""
    
    content = snakefile_path.read_text()
    rules_content = rules_path.read_text()
    
    # Find the insertion point (after generate_consensus rule)
    pattern = r'(rule generate_consensus:.*?""")\n\n(# Rule 8:)'
    
    replacement = r'\1\n\n' + rules_content + r'\n\n\2'
    
    new_content = re.sub(pattern, replacement, content, flags=re.DOTALL, count=1)
    
    if new_content == content:
        print("❌ Could not find insertion point in Snakefile", file=sys.stderr)
        return False
    
    snakefile_path.write_text(new_content)
    print(f"✅ Inserted polishing rules")
    return True


def update_downstream_references(snakefile_path: Path) -> int:
    """Update downstream rules to use polished consensus."""
    
    content = snakefile_path.read_text()
    
    # Pattern to replace: rules.generate_consensus.output.fasta
    # With: rules.polish_consensus.output.polished
    pattern = r'rules\.generate_consensus\.output\.fasta'
    replacement = r'rules.polish_consensus.output.polished'
    
    new_content = re.sub(pattern, replacement, content)
    count = len(re.findall(pattern, content))
    
    snakefile_path.write_text(new_content)
    print(f"✅ Updated {count} downstream references")
    return count


def update_blast_novelty_input(snakefile_path: Path) -> bool:
    """Update blast_novelty_scan to use polished consensus."""
    
    content = snakefile_path.read_text()
    
    # Find and replace the blast_novelty_scan consensus input
    pattern = r'(rule blast_novelty_scan:.*?input:.*?consensus=)"\{output_dir\}/\{\{sample\}\}/09_consensus/\{\{sample\}\}_consensus\.fasta"\.format\([^)]+\)'
    replacement = r'\1rules.polish_consensus.output.polished'
    
    new_content = re.sub(pattern, replacement, content, flags=re.DOTALL, count=1)
    
    if new_content != content:
        snakefile_path.write_text(new_content)
        print(f"✅ Updated blast_novelty_scan input")
        return True
    else:
        print("⚠️  Could not update blast_novelty_scan (may already be updated)")
        return False


def add_config_parameter(config_path: Path) -> bool:
    """Add pipeline_mode to config.yaml."""
    
    content = config_path.read_text()
    
    # Check if already exists
    if 'pipeline_mode' in content:
        print("ℹ️  pipeline_mode already in config")
        return True
    
    # Add after the Tier 3 section
    addition = '\n# Consensus generation mode\npipeline_mode: "LABORATORY_FULL"  # or "FIELD_RAPID"\n'
    
    # Find Tier 3 section
    if 'Tier 3' in content:
        pattern = r'(# Tier 3.*?\n.*?\n)'
        new_content = re.sub(pattern, r'\1' + addition, content, flags=re.DOTALL)
        config_path.write_text(new_content)
        print("✅ Added pipeline_mode to config.yaml")
        return True
    else:
        # Append to end
        config_path.write_text(content + addition)
        print("✅ Added pipeline_mode to config.yaml (at end)")
        return True


def main():
    print("=" * 60)
    print("🔧 Vibrion Pipeline: Polishing Integration")
    print("=" * 60)
    
    # Paths
    base_dir = Path(__file__).parent.parent
    snakefile = base_dir / "workflow" / "Snakefile"
    rules_file = base_dir / "workflow" / "rules" / "polishing.smk"
    config_file = base_dir / "workflow" / "config" / "config.yaml"
    
    # Verify files exist
    if not snakefile.exists():
        print(f"❌ Snakefile not found: {snakefile}", file=sys.stderr)
        sys.exit(1)
    
    if not rules_file.exists():
        print(f"❌ Polishing rules not found: {rules_file}", file=sys.stderr)
        sys.exit(1)
    
    # Create backup
    print("\n📁 Creating backups...")
    backup_file(snakefile)
    backup_file(config_file)
    
    # Step 1: Insert polishing rules
    print("\n🔧 Step 1: Inserting polishing rules...")
    if not insert_polishing_rules(snakefile, rules_file):
        print("❌ Integration failed", file=sys.stderr)
        sys.exit(1)
    
    # Step 2: Update downstream references
    print("\n🔧 Step 2: Updating downstream references...")
    count = update_downstream_references(snakefile)
    
    # Step 3: Update BLAST input
    print("\n🔧 Step 3: Updating BLAST novelty scan...")
    update_blast_novelty_input(snakefile)
    
    # Step 4: Update config
    print("\n🔧 Step 4: Updating config.yaml...")
    add_config_parameter(config_file)
    
    print("\n" + "=" * 60)
    print("✅ Integration Complete!")
    print("=" * 60)
    print("\n📋 Next steps:")
    print("1. Review changes: git diff workflow/Snakefile")
    print("2. Test dry-run: snakemake -n --snakefile workflow/Snakefile")
    print("3. If issues occur, restore from .backup files")
    print("\n🎯 New rules added: detect_platform, polish_consensus")


if __name__ == "__main__":
    main()
