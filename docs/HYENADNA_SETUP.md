# HyenaDNA Installation Guide

HyenaDNA is used in Vibrion Sentinel's **Tier 1 local triage** step to detect genomic anomalies without cloud calls. This guide covers installation on macOS (Apple Silicon + Intel) and Linux.

---

## What it's used for

The `local_triage` rule loads HyenaDNA to compute per-locus anomaly scores against
the Haiti 2010 and 2022 baselines. If the score is below the threshold the sample is
classified **LOCAL VERIFIED SAFE** and cloud Evo2 escalation is skipped entirely,
saving API quota and reducing latency.

Model used: `LongSafari/hyenadna-tiny-1k-seqlen-hf` (~170MB, runs on CPU)

---

## Requirements

| Component | Version |
|-----------|---------|
| Python | ≥ 3.9 |
| PyTorch | ≥ 2.0 |
| Transformers (HuggingFace) | ≥ 4.35 |
| einops | ≥ 0.6 |
| Flash-Attn | optional (GPU only) |

---

## Installation

### Option A — via the vibrion conda environment (recommended)

```bash
conda env create -f environment.yml
conda activate vibrion
```

The `environment.yml` includes all HyenaDNA dependencies. On first run the pipeline
will automatically download the model weights to `data/models/triage/`.

### Option B — manual pip install

```bash
pip install transformers>=4.35 einops torch torchvision torchaudio
```

### Option C — Apple Silicon (M1/M2/M3) with MPS acceleration

```bash
pip install torch torchvision torchaudio  # ships with MPS backend
pip install transformers einops
```

Set in your config YAML:
```yaml
triage:
  device: "mps"   # or "cpu" for Intel Mac / Linux without GPU
```

---

## Verify installation

```python
from transformers import AutoTokenizer, AutoModelForSequenceClassification
model = AutoTokenizer.from_pretrained("LongSafari/hyenadna-tiny-1k-seqlen-hf",
                                       trust_remote_code=True)
print("HyenaDNA tokenizer loaded ✅")
```

---

## Downloading the model weights manually (offline environments)

```bash
pip install huggingface_hub
python3 -c "
from huggingface_hub import snapshot_download
snapshot_download(
    repo_id='LongSafari/hyenadna-tiny-1k-seqlen-hf',
    local_dir='data/models/triage/hyenadna-tiny-1k-seqlen-hf',
    ignore_patterns=['*.gguf']
)
print('Done')
"
```

Point the config to the local path:
```yaml
triage:
  hyena_model: "data/models/triage/hyenadna-tiny-1k-seqlen-hf"
```

---

## Troubleshooting

| Error | Fix |
|-------|-----|
| `ModuleNotFoundError: einops` | `pip install einops` |
| `trust_remote_code` warning | Pass `trust_remote_code=True` to `from_pretrained` — expected for HyenaDNA |
| `MPS backend out of memory` | Set `device: "cpu"` in triage config |
| Model not found at runtime | Check `triage.hyena_model` in your config YAML |
| Slow on CPU | Normal — tiny model takes ~2-5s per locus on CPU. Use `hyena_use_real_model: false` to skip and use k-mer fallback |

---

## Disabling HyenaDNA (k-mer fallback mode)

If you cannot install PyTorch or are on a restricted system:

```yaml
triage:
  hyena_use_real_model: false   # Falls back to sourmash k-mer comparison only
```

The pipeline will still run Tier 0 (sourmash) triage. Only Tier 1 anomaly scoring
is skipped. Results are slightly less sensitive to novel variants.
