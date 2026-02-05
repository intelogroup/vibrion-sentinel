import sys
from pathlib import Path
import json
import torch
import torch.nn as nn
import torch.nn.functional as F
import argparse
import numpy as np
from Bio import SeqIO
import math

# Add HyenaDNA repo to path
# Add HyenaDNA repo to path
REPO_PATH = Path("data/external/hyena-dna").resolve()
sys.path.append(str(REPO_PATH))

# Import Hyena modules (HyenaOperator is safe from flash_attn)
try:
    from src.models.sequence.hyena import HyenaOperator
except ImportError as e:
    print(f"❌ Could not import HyenaDNA modules: {e}")
    exit(1)

# Simple Tokenizer to bypass broken repo version
class SimpleCharTokenizer:
    def __init__(self, characters=['A', 'C', 'G', 'T', 'N']):
        self.char_to_id = {ch: i + 7 for i, ch in enumerate(characters)}
        self.char_to_id["[UNK]"] = 6
        self.unk_id = 6

    def encode(self, text):
        return [self.char_to_id.get(ch, self.unk_id) for ch in text]

# Device Configuration
device = "mps" if torch.backends.mps.is_available() else "cpu"
print(f"🦁 Device: {device}")

# --- Standalone Model Definitions [Bypassing flash_attn] ---

class StandaloneMlp(nn.Module):
    def __init__(self, d_model, d_inner, activation_fn=F.gelu):
        super().__init__()
        self.fc1 = nn.Linear(d_model, d_inner)
        self.fc2 = nn.Linear(d_inner, d_model)
        self.act = activation_fn
        
    def forward(self, x):
        return self.fc2(self.act(self.fc1(x)))

class StandaloneBlock(nn.Module):
    def __init__(self, config, layer_idx):
        super().__init__()
        d_model = config['d_model']
        d_inner = config['d_inner']
        
        # Mixer (Hyena)
        # Config has "layer": {"_name_": "hyena", ...}
        hyena_cfg = config['layer']
        # Extract args for HyenaOperator
        self.mixer = HyenaOperator(
            d_model=d_model,
            l_max=hyena_cfg['l_max'],
            order=hyena_cfg.get('order', 2),
            filter_order=hyena_cfg.get('filter_order', 64),
            num_heads=hyena_cfg.get('num_heads', 1), # Default 1
            inner_factor=hyena_cfg.get('inner_factor', 1),
            modulate=hyena_cfg.get('modulate', True),
            w=hyena_cfg.get('w', 10), # from checked config
            lr=hyena_cfg.get('lr', 1e-3),
            lr_pos_emb=hyena_cfg.get('lr_pos_emb', 1e-5),
            dropout=hyena_cfg.get('dropout', 0.0),
            emb_dim=hyena_cfg.get('emb_dim', 3) # Fix: Pass emb_dim (default 3, ckpt uses 5)
        )

        
        # MLP
        self.mlp = StandaloneMlp(d_model, d_inner)
        
        # Layer Norms
        self.norm1 = nn.LayerNorm(d_model, eps=1e-5)
        self.norm2 = nn.LayerNorm(d_model, eps=1e-5)
        
        # Dropout (simplified)
        self.resid_dropout = nn.Dropout(config.get('resid_dropout', 0.0))

    def forward(self, hidden_states):
        # Hyena Path
        residual = hidden_states
        hidden_states = self.norm1(hidden_states)
        
        # Mixer expect (B, L, D) -> (B, D, L) usually? 
        # Check HyenaOperator.forward: u = self.in_proj(u) -> rearrange b l d -> b d l
        # So it expects (B, L, D).
        mx_out = self.mixer(hidden_states)
        
        hidden_states = residual + self.resid_dropout(mx_out)
        
        # MLP Path
        residual = hidden_states
        hidden_states = self.norm2(hidden_states)
        mlp_out = self.mlp(hidden_states)
        
        hidden_states = residual + self.resid_dropout(mlp_out)
        
        return hidden_states

class StandaloneBackbone(nn.Module):
    def __init__(self, config):
        super().__init__()
        d_model = config['d_model']
        vocab_size = config['vocab_size']
        n_layer = config['n_layer']
        
        # Embeddings
        # We assume no positional embeddings based on checks
        self.embeddings = nn.Module()
        self.embeddings.word_embeddings = nn.Embedding(vocab_size, d_model)
        
        # Layers
        self.layers = nn.ModuleList([
            StandaloneBlock(config, i) for i in range(n_layer)
        ])
        
        # Final Norm
        self.ln_f = nn.LayerNorm(d_model, eps=1e-5)
        
    def forward(self, input_ids):
        # Embed
        hidden_states = self.embeddings.word_embeddings(input_ids)
        
        # Pass through layers
        for layer in self.layers:
            hidden_states = layer(hidden_states)
            
        return self.ln_f(hidden_states)

class StandaloneHyenaLM(nn.Module):
    def __init__(self, config):
        super().__init__()
        d_model = config['d_model']
        
        # Handle Vocab Padding (Config says 12, Ckpt has 16 due to pad_vocab_size_multiple=8)
        vocab_size = config['vocab_size']
        if config.get('pad_vocab_size_multiple', 1) > 1:
            vocab_size = (math.ceil(vocab_size / config['pad_vocab_size_multiple']) * config['pad_vocab_size_multiple'])
            
        self.backbone = StandaloneBackbone(config)
        # Update backbone embeddings to use padded vocab size
        self.backbone.embeddings.word_embeddings = nn.Embedding(vocab_size, d_model)
        
        self.lm_head = nn.Linear(d_model, vocab_size, bias=False)
        
        # Tie weights
        self.lm_head.weight = self.backbone.embeddings.word_embeddings.weight

    def forward(self, input_ids, return_embeddings=False):
        hidden_states = self.backbone(input_ids)
        if return_embeddings:
            return hidden_states
        logits = self.lm_head(hidden_states)
        return logits

def load_hyena_model_standalone(model_dir):
    """
    Load weights into StandaloneHyenaLM
    """
    model_dir = Path(model_dir)
    config_path = model_dir / "config.json"
    ckpt_path = model_dir / "weights.ckpt"
    
    with open(config_path) as f:
        config = json.load(f)
        
    print("   Instantiating StandaloneHyenaLM...")
    model = StandaloneHyenaLM(config)
    
    print(f"   Loading Weights: {ckpt_path}")
    state_dict = torch.load(ckpt_path, map_location='cpu', weights_only=False)
    if "state_dict" in state_dict:
        state_dict = state_dict["state_dict"]
        
    # Remap keys
    new_state_dict = {}
    for k, v in state_dict.items():
        # Strip 'model.' prefix
        if k.startswith("model."):
            k = k[6:]
            
        # Strip '.layer' from mixer/mlp if present (due to CheckpointedModule)
        # e.g. backbone.layers.0.mixer.layer.filter_fn -> backbone.layers.0.mixer.filter_fn
        k = k.replace("mixer.layer.", "mixer.")
        k = k.replace("mlp.layer.", "mlp.")
        
        # Map embedding keys if valid
        # backbone.embeddings.word_embeddings.weight matches
        
        # Filter torchmetrics
        if "torchmetrics" in k:
            continue
            
        new_state_dict[k] = v
        
    msg = model.load_state_dict(new_state_dict, strict=False)
    print(f"   Weights Loaded. Missing: {len(msg.missing_keys)}, Unexpected: {len(msg.unexpected_keys)}")
    if len(msg.missing_keys) > 0:
        print(f"   Sample Missing: {msg.missing_keys[:5]}")
    
    model.to(device)
    model.eval()
    return model, config

def get_embedding(model, tokenizer, sequence):
    tokens = tokenizer.encode(sequence)
    input_ids = torch.tensor([tokens]).to(device)
    with torch.no_grad():
        # Get hidden states (B, L, D)
        embeddings = model(input_ids, return_embeddings=True)
        # Mean pooling (B, D)
        mean_embedding = embeddings.mean(dim=1).squeeze().cpu().numpy()
    
    del input_ids, embeddings
    if device == "mps":
        torch.mps.empty_cache()
        
    return mean_embedding

def main():
    parser = argparse.ArgumentParser(description="Run HyenaDNA 1M Analysis")
    parser.add_argument("--input", required=True, help="Input FASTA")
    parser.add_argument("--model-dir", required=True, help="Path to model weights")
    parser.add_argument("--output", required=True, help="Output JSON")
    parser.add_argument("--reference", help="Optional Reference FASTA for Cosine Similarity")
    args = parser.parse_args()
    
    # 1. Load Model
    model, config = load_hyena_model_standalone(args.model_dir)
    
    # 2. Tokenizer (Simple standalone)
    tokenizer = SimpleCharTokenizer()
    
    # 2b. Pre-compute Reference Embeddings if provided
    ref_embeddings = {}
    if args.reference:
        print(f"   Loading Reference: {args.reference}")
        ref_records = list(SeqIO.parse(args.reference, "fasta"))
        for rec in ref_records:
            ref_embeddings[rec.id] = get_embedding(model, tokenizer, str(rec.seq).upper())
        print(f"   Computed {len(ref_embeddings)} reference embeddings.")

    # 3. Read Sequences
    records = list(SeqIO.parse(args.input, "fasta"))
    print(f"   Input FASTA contains {len(records)} records.")
    
    all_results = {}
    
    # 4. Inference loop
    for record in records:
        sequence = str(record.seq).upper()
        seq_len = len(sequence)
        locus_name = record.id
        print(f"   Processing Locus: {locus_name} ({seq_len} bp)")
        
        # Inference (Logits for PPL)
        tokens = tokenizer.encode(sequence)
        input_ids = torch.tensor([tokens]).to(device)
        
        with torch.no_grad():
            logits = model(input_ids)
            
        # PPL
        shift_logits = logits[..., :-1, :].contiguous()
        shift_labels = input_ids[..., 1:].contiguous()
        loss_fct = torch.nn.CrossEntropyLoss(reduction='none')
        loss = loss_fct(shift_logits.view(-1, shift_logits.size(-1)), shift_labels.view(-1))
        
        avg_loss = loss.mean().item()
        perplexity = np.exp(avg_loss)
        
        # Similarity
        cosine_sim = None
        if locus_name in ref_embeddings:
            # Get input embedding
            # Note: We computed logits above, but need hidden states now.
            # To save memory/compute, we could have done both, but forward logic is split.
            # Rerunning for embedding is safer for memory than holding both.
            # Actually, let's just run embedding generation.
            
            # Clean up logits first
            del logits, input_ids
            if device == "mps": torch.mps.empty_cache()
            
            input_emb = get_embedding(model, tokenizer, sequence)
            ref_emb = ref_embeddings[locus_name]
            
            # Cosine Similarity
            dot_product = np.dot(input_emb, ref_emb)
            norm_a = np.linalg.norm(input_emb)
            norm_b = np.linalg.norm(ref_emb)
            cosine_sim = float(dot_product / (norm_a * norm_b))
            print(f"      Loss={avg_loss:.4f} PPL={perplexity:.2f} Sim={cosine_sim:.4f}")
        else:
            print(f"      Loss={avg_loss:.4f} PPL={perplexity:.2f}")
            # Clean up logits
            del logits, input_ids
            if device == "mps": torch.mps.empty_cache()
        
        all_results[locus_name] = {
            "length": seq_len,
            "loss": avg_loss,
            "perplexity": perplexity,
            "cosine_similarity": cosine_sim
        }
            
    # 5. Save
    with open(args.output, "w") as f:
        json.dump({
            "sample": Path(args.input).stem,
            "model": "HyenaDNA-1M",
            "device": device,
            "loci_analysis": all_results
        }, f, indent=2)
        
    print(f"✅ Analysis saved to {args.output}")

if __name__ == "__main__":
    main()
