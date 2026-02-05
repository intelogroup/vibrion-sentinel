#!/usr/bin/env python3
"""
Vibrion Sentinel Memory Manager
Non-blocking adaptive memory management system for genomic analysis pipelines.

This module provides transparent, non-intrusive memory monitoring and adaptive
resource allocation without interrupting pipeline execution. All decisions are
informational and reversible.

Key Design Principles:
1. NON-BLOCKING: Never stops pipeline, only adapts configuration
2. SILENT: Logs decisions but doesn't interrupt with warnings
3. GRACEFUL: Drops features (phylogeny, visualization) before core analysis
4. OBSERVABLE: All decisions logged and reviewable
5. REVERSIBLE: No permanent changes, only config adaptation
"""

import json
import psutil
import logging
import threading
import time
from typing import Dict, Optional
from pathlib import Path
from datetime import datetime


class MemoryManager:
    """Adaptive memory management for Snakemake workflows."""

    # Memory tier thresholds (GB)
    TIER_THRESHOLDS = {
        'FULL': 10,
        'BALANCED': 6,
        'BUNKER': 3,
        'EMERGENCY': 0
    }

    # Feature availability by tier (lower tiers include features from higher tiers)
    TIER_FEATURES = {
        'FULL': {
            'kraken_db': 'data/kraken2_standard_8gb',
            'threads': 4,
            'spades_kmers': '21,33,55,77',
            'mafft_strategy': 'full',
            'mafft_max_refs': None,
            'features': ['kraken2_full', 'mmseqs2_rescue', 'mafft_full', 'evo2', 'phylogeny']
        },
        'BALANCED': {
            'kraken_db': 'data/kraken2_haiti_custom',
            'threads': 3,
            'spades_kmers': '21,33,55',
            'mafft_strategy': 'fast',
            'mafft_max_refs': 50,
            'features': ['kraken2_haiti', 'mmseqs2_rescue', 'mafft_fast', 'evo2']
        },
        'BUNKER': {
            'kraken_db': 'data/kraken2_haiti_custom',
            'threads': 2,
            'spades_kmers': '21,33',
            'mafft_strategy': 'skip',
            'mafft_max_refs': 0,
            'features': ['kraken2_haiti', 'mmseqs2_rescue']
        },
        'EMERGENCY': {
            'kraken_db': 'data/kraken2_haiti_custom',
            'threads': 1,
            'spades_assembly': False,
            'mafft_strategy': 'skip',
            'features': ['kraken2_haiti']
        }
    }

    def __init__(self, config: Dict, enable_monitoring: bool = True):
        """
        Initialize memory manager.

        Args:
            config: Snakemake config dict (passed via configfile)
            enable_monitoring: Enable background memory monitoring
        """
        self.config = config
        self.memory_config = config.get('memory_management', {})
        self.current_tier = None
        self.monitoring_enabled = enable_monitoring
        self.monitor_thread = None
        self.monitoring_stop = threading.Event()

        # Setup logging
        self.logger = self._setup_logging()
        self.logger.info("MemoryManager initialized")

        # Determine initial tier
        self.current_tier = self._detect_tier()
        available_gb = self._get_available_memory_gb()
        self.logger.info(f"Memory tier: {self.current_tier} "
                         f"(available: {available_gb:.1f}GB)")

        # Start background monitoring if enabled
        if self.memory_config.get('monitoring', {}).get('enabled', False):
            self._start_monitoring()

    def _setup_logging(self) -> logging.Logger:
        """Configure logging for memory manager."""
        logger = logging.getLogger('vibrion.memory')
        logger.setLevel(logging.INFO)

        log_dir = Path('logs')
        log_dir.mkdir(exist_ok=True)

        memory_monitoring_config = self.memory_config.get('monitoring', {})
        log_file_name = memory_monitoring_config.get('log_file', 'memory_profile.log')

        # Handle nested paths (e.g., 'logs/memory_profile.log')
        log_file = Path(log_file_name)
        if log_file.parent != Path('.'):
            # If specified path has parent dir, create it
            log_file.parent.mkdir(parents=True, exist_ok=True)
        else:
            # Otherwise put in logs/ directory
            log_file = log_dir / log_file_name

        handler = logging.FileHandler(log_file)
        formatter = logging.Formatter(
            '%(asctime)s [%(levelname)s] %(name)s: %(message)s'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        return logger

    def _get_available_memory_gb(self) -> float:
        """Get available system memory in GB."""
        return psutil.virtual_memory().available / (1024 ** 3)

    def _detect_tier(self, override_tier: Optional[str] = None) -> str:
        """
        Detect appropriate memory tier.

        Args:
            override_tier: Force specific tier (for testing)

        Returns:
            Tier name (FULL, BALANCED, BUNKER, EMERGENCY)
        """
        # Check for forced tier in config
        if self.memory_config.get('force_tier'):
            tier = self.memory_config['force_tier']
            self.logger.info(f"Using forced tier: {tier}")
            return tier

        # Auto-detect based on available memory
        available_gb = self._get_available_memory_gb()

        for tier in ['FULL', 'BALANCED', 'BUNKER', 'EMERGENCY']:
            if available_gb >= self.TIER_THRESHOLDS[tier]:
                return tier

        return 'EMERGENCY'  # Fallback

    def get_config_for_tier(self, tier: Optional[str] = None) -> Dict:
        """
        Get configuration parameters for a specific tier.

        Args:
            tier: Tier name (uses current if not specified)

        Returns:
            Dict of config overrides for this tier
        """
        tier = tier or self.current_tier
        return self.TIER_FEATURES.get(tier, {})

    def get_effective_config(self) -> Dict:
        """Get effective Snakemake config with memory-aware overrides."""
        tier_config = self.get_config_for_tier()

        # Build effective config
        effective = {
            'threads': tier_config.get('threads', self.config.get('threads', 4)),
            'kraken_db': tier_config.get('kraken_db', self.config.get('kraken_db')),
            'spades_kmers': tier_config.get('spades_kmers'),
            'mafft_strategy': tier_config.get('mafft_strategy', 'full'),
            'mafft_max_refs': tier_config.get('mafft_max_refs'),
            'memory_tier': self.current_tier,
            'enabled_features': tier_config.get('features', [])
        }

        return effective

    def should_run_feature(self, feature_name: str) -> bool:
        """Check if feature is enabled in current tier."""
        tier_config = self.get_config_for_tier()
        return feature_name in tier_config.get('features', [])

    def _start_monitoring(self):
        """Start background memory monitoring thread."""
        self.monitor_thread = threading.Thread(
            target=self._monitor_memory,
            daemon=True
        )
        self.monitor_thread.start()
        self.logger.info("Background memory monitoring started")

    def _monitor_memory(self):
        """
        Background monitoring thread (non-blocking).

        Checks memory usage periodically and logs alerts if threshold exceeded.
        Never stops pipeline, only generates informational logs.
        """
        alert_threshold = self.memory_config.get('monitoring', {}).get(
            'alert_threshold_percent', 85
        )
        check_interval = self.memory_config.get('monitoring', {}).get(
            'check_interval_sec', 30
        )

        while not self.monitoring_stop.is_set():
            try:
                memory_info = psutil.virtual_memory()
                usage_percent = memory_info.percent

                # Check for tier change (non-blocking adaptation)
                previous_tier = self.current_tier
                available_gb = memory_info.available / (1024 ** 3)

                # Only log if significant change
                if usage_percent > alert_threshold:
                    used_gb = memory_info.used / (1024**3)
                    total_gb = memory_info.total / (1024**3)
                    self.logger.warning(
                        f"Memory usage high: {usage_percent:.1f}% "
                        f"({used_gb:.1f}GB / {total_gb:.1f}GB)"
                    )

                # Check if tier has changed
                new_tier = self._detect_tier()
                if new_tier != previous_tier:
                    self.current_tier = new_tier
                    self.logger.info(
                        f"Tier change detected: {previous_tier} → {new_tier} "
                        f"(available: {available_gb:.1f}GB)"
                    )

                time.sleep(check_interval)
            except Exception as e:
                self.logger.error(f"Monitoring error: {e}")
                time.sleep(check_interval)

    def stop_monitoring(self):
        """Stop background monitoring thread."""
        if self.monitor_thread:
            self.monitoring_stop.set()
            self.monitor_thread.join(timeout=5)
            self.logger.info("Background memory monitoring stopped")

    def get_memory_summary(self) -> Dict:
        """Get current memory state summary."""
        memory = psutil.virtual_memory()
        swap = psutil.swap_memory()

        return {
            'timestamp': datetime.now().isoformat(),
            'tier': self.current_tier,
            'available_gb': memory.available / (1024 ** 3),
            'used_gb': memory.used / (1024 ** 3),
            'total_gb': memory.total / (1024 ** 3),
            'percent': memory.percent,
            'swap_used_gb': swap.used / (1024 ** 3),
            'swap_total_gb': swap.total / (1024 ** 3),
            'enabled_features': self.TIER_FEATURES[self.current_tier].get('features', [])
        }

    def write_memory_profile(self, sample_id: str, output_dir: str = 'data/pipeline_output'):
        """Write memory profile to JSON file."""
        profile = self.get_memory_summary()
        profile['sample_id'] = sample_id

        output_path = Path(output_dir) / sample_id / '00_metadata'
        output_path.mkdir(parents=True, exist_ok=True)

        profile_file = output_path / 'memory_profile.json'
        with open(profile_file, 'w') as f:
            json.dump(profile, f, indent=2)

        self.logger.info(f"Memory profile written to {profile_file}")


# Global instance (initialized by Snakefile)
_memory_manager: Optional[MemoryManager] = None


def initialize_memory_manager(snakemake_config: Dict) -> MemoryManager:
    """
    Initialize global memory manager (called from Snakefile).

    Args:
        snakemake_config: Snakemake config dict

    Returns:
        MemoryManager instance
    """
    global _memory_manager
    if _memory_manager is None:
        _memory_manager = MemoryManager(snakemake_config)
    return _memory_manager


def get_memory_manager() -> Optional[MemoryManager]:
    """Get global memory manager instance."""
    return _memory_manager


def cleanup_memory_manager():
    """Clean up memory manager on exit."""
    global _memory_manager
    if _memory_manager:
        _memory_manager.stop_monitoring()
        _memory_manager = None


if __name__ == '__main__':
    # Standalone testing
    import yaml

    # Load config for testing
    config_path = 'workflow/config/config.yaml'
    with open(config_path) as f:
        test_config = yaml.safe_load(f)

    # Initialize and show summary
    mm = MemoryManager(test_config, enable_monitoring=True)

    print("\n=== Memory Manager Status ===")
    summary = mm.get_memory_summary()
    for key, value in summary.items():
        print(f"{key}: {value}")

    print("\nEffective Configuration:")
    effective = mm.get_effective_config()
    for key, value in effective.items():
        print(f"  {key}: {value}")

    # Clean up
    mm.stop_monitoring()
