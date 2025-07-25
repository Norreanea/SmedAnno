#!/usr/bin/env python3
"""
DeepSplice (SmedAnno fork) — CUDA-ready, batched
--------------------------------------------------------------
"""
from __future__ import annotations

import argparse, bz2, gzip, sys, warnings
from pathlib import Path
from typing import Dict, Iterator, Tuple, List

import numpy as np
import torch
import torch.nn as nn
from Bio import SeqIO

# ────────────────────────────── Constants ──────────────────────────────
WINDOW = 10_240
STEP   = WINDOW // 2
DEFAULT_BATCH = 16  # safe for 2 GB GT 1030; tune upward if you have more VRAM

CUDA  = torch.cuda.is_available()
DEVICE = "cuda" if CUDA else "cpu"
if CUDA:
    torch.backends.cudnn.benchmark = True  # let cuDNN pick fastest conv algos

SPECIES = {"human", "mouse", "zebrafish", "honeybee", "thalecress"}

# ────────────────────────────── Model blocks ────────────────────────────
class ResidualUnit(nn.Module):
    def __init__(self, channels: int, kernel_size: int):
        super().__init__()
        pad = kernel_size // 2
        self.bn1 = nn.BatchNorm1d(channels)
        self.act = nn.GELU()
        self.conv1 = nn.Conv1d(channels, channels, kernel_size, padding=pad)
        self.bn2 = nn.BatchNorm1d(channels)
        self.conv2 = nn.Conv1d(channels, channels, kernel_size, padding=pad)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        y = self.conv1(self.act(self.bn1(x)))
        y = self.conv2(self.act(self.bn2(y)))
        return x + y

class DilatedConvolutionalUnit(nn.Module):
    def __init__(self, channels: int, kernel_size: int, dilation: int):
        super().__init__()
        pad = (kernel_size - 1) * dilation // 2
        self.conv = nn.Conv1d(channels, channels, kernel_size, padding=pad, dilation=dilation)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.conv(x)

class OSAI(nn.Module):
    """OpenSpliceAI-like architecture used by DeepSplice.

    IMPORTANT:
    - `initial_conv` consumes 4 channels (one-hot DNA) → CHANNELS
    - `initial_skip` consumes the already projected CHANNELS tensor
    This matches old checkpoints that had 32→32 weights for skip.
    """
    def __init__(self):
        super().__init__()
        CHANNELS = 32
        # project one-hot to feature space
        self.initial_conv = nn.Sequential(nn.Conv1d(4, CHANNELS, kernel_size=1))
        # skip path operates on CHANNELS already
        self.initial_skip = nn.Sequential(nn.Conv1d(CHANNELS, CHANNELS, kernel_size=1))

        units: List[nn.Module] = []
        for i in range(20):
            if i in {0,1,2,3,5,6,7,8}:
                units.append(ResidualUnit(CHANNELS, 11))
            elif i in {4,9,14,19}:  # dilated
                dilation = 2 ** (i // 5 + 1) if i != 19 else 25
                units.append(DilatedConvolutionalUnit(CHANNELS, 1, dilation))
            elif i in {10,11,12,13}:
                units.append(ResidualUnit(CHANNELS, 21))
            elif i in {15,16,17,18}:
                units.append(ResidualUnit(CHANNELS, 41))
        self.residual_units = nn.Sequential(*units)
        self.final_conv = nn.Conv1d(CHANNELS, 3, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x0 = self.initial_conv(x)
        x  = self.residual_units(x0 + self.initial_skip(x0))
        return torch.softmax(self.final_conv(x), dim=1)

# ────────────────────────────── I/O helpers ─────────────────────────────

def one_hot_encode(seq_str: str) -> np.ndarray:
    seq = np.array(list(seq_str.upper()))
    mapping = np.array(['A','C','G','T'])
    one_hot = (seq[..., None] == mapping).astype(np.uint8)
    one_hot[seq == 'N'] = 0
    return one_hot.T  # (4, L)

def fasta_iter(handle) -> Iterator[Tuple[str, np.ndarray]]:
    for rec in SeqIO.parse(handle, "fasta"):
        yield rec.id, one_hot_encode(str(rec.seq))

def sliding_windows(seq: np.ndarray):
    L = seq.shape[1]
    for off in range(0, L, STEP):
        chunk = seq[:, off:off+WINDOW]
        if chunk.shape[1] < WINDOW:
            pad = WINDOW - chunk.shape[1]
            chunk = np.pad(chunk, ((0,0),(0,pad)), "constant")
        yield off, chunk

# ────────────────────────────── Peak caller ─────────────────────────────

def call_peaks(track: torch.Tensor, chrom: str, offset: int, kind: str, thr: float, out) -> None:
    mask = track > thr
    prev = torch.zeros_like(mask)
    prev[1:] = mask[:-1]
    starts = torch.nonzero(mask & ~prev, as_tuple=True)[0]
    for pos in starts.tolist():
        prob = float(track[pos])
        coord = offset + pos + 1  # GFF3 is 1-based
        out.write(f"{chrom}\tDeepSplice\t{kind}\t{coord}\t{coord}\t{prob:.3f}\t.\t.\t.\n")

# ────────────────────────────── Checkpoint loader ───────────────────────

def _safe_torch_load(path: Path, device: str) -> Dict[str, torch.Tensor]:
    try:
        return torch.load(path, map_location=device, weights_only=True)
    except TypeError:
        warnings.warn("torch.load(weights_only=True) unsupported; falling back to pickle.")
        return torch.load(path, map_location=device)

def _extract_state_dict(raw):
    if isinstance(raw, dict):
        for k in ("state_dict","model_state_dict","model","net"):
            if k in raw and isinstance(raw[k], dict):
                return raw[k]
    return raw

def _rename_legacy_keys(state: Dict[str, torch.Tensor], target_keys) -> Dict[str, torch.Tensor]:
    has_seq_skip = any(k.startswith("initial_skip.0.") for k in target_keys)
    has_seq_conv = any(k.startswith("initial_conv.0.") for k in target_keys)
    out = {}
    for k,v in state.items():
        nk = k
        for pref in ("module.","model."):
            if nk.startswith(pref):
                nk = nk[len(pref):]
        nk = nk.replace(".conv.", ".")
        if nk.startswith("initial_skip.") and has_seq_skip and not nk.startswith("initial_skip.0."):
            nk = nk.replace("initial_skip.", "initial_skip.0.")
        if nk.startswith("initial_conv.") and has_seq_conv and not nk.startswith("initial_conv.0."):
            nk = nk.replace("initial_conv.", "initial_conv.0.")
        out[nk] = v
    return out

def load_checkpoint_into(model: nn.Module, ckpt_path: Path, device: str) -> None:
    raw   = _safe_torch_load(ckpt_path, device)
    state = _extract_state_dict(raw)
    if not isinstance(state, dict):
        raise RuntimeError(f"Unsupported checkpoint format: {ckpt_path}")

    fixed = _rename_legacy_keys(state, model.state_dict().keys())
    # Drop any tensor whose shape doesn't match
    model_sd = model.state_dict()
    good, skipped = {}, []
    for k,v in fixed.items():
        if k in model_sd and model_sd[k].shape == v.shape:
            good[k] = v
        else:
            skipped.append(k)
    missing, unexpected = model.load_state_dict(good, strict=False)
    if skipped or missing or unexpected:
        print(f"[DeepSplice] Warning: skipped(shape)={skipped}, missing={list(missing)}, unexpected={list(unexpected)}", file=sys.stderr)

# ────────────────────────────── Batch runner ────────────────────────────

def _run_batch(model: nn.Module, chunks: List[torch.Tensor], offsets: List[int], chrom: str, thr: float, gff):
    batch = torch.stack(chunks).to(DEVICE, non_blocking=True)
    with torch.inference_mode():
        if DEVICE == "cuda":
            with torch.autocast("cuda", dtype=torch.float16):
                probs = model(batch).float().cpu()
        else:
            probs = model(batch).cpu()
    for i, off in enumerate(offsets):
        call_peaks(probs[i, 1], chrom, off, "donor",    thr, gff)
        call_peaks(probs[i, 2], chrom, off, "acceptor", thr, gff)

# ────────────────────────────── CLI ─────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser("OpenSpliceAI splice-site predictor")
    p.add_argument("fasta", help="Genome FASTA (optionally gz/bz2)")
    p.add_argument("-o", "--out", default="splice_predictions.gff3", help="Output GFF3")
    p.add_argument("--thr", type=float, default=0.95, help="Posterior threshold")
    p.add_argument("--species", choices=SPECIES, help="Shortcut to DSmodels/model_<species>.pt")
    p.add_argument("--weights", help="Path to custom *.pt checkpoint")
    p.add_argument("--id-prefix", default="", help="Prefix for seq IDs in GFF3")
    p.add_argument("--batch", type=int, default=DEFAULT_BATCH, help="Batch size for CUDA/CPU")
    return p.parse_args()

# ────────────────────────────── Main ────────────────────────────────────

def main():
    args = parse_args()

    # choose checkpoint
    if args.weights:
        ckpt = Path(args.weights)
    elif args.species:
        ckpt = Path(__file__).resolve().parent / "DSmodels" / f"model_{args.species}.pt"
    else:
        sys.exit("[DeepSplice] Provide --species or --weights")
    if not ckpt.exists():
        sys.exit(f"[DeepSplice] Checkpoint not found: {ckpt}")

    model = OSAI().to(DEVICE)
    load_checkpoint_into(model, ckpt, DEVICE)
    model.eval()

    opener = gzip.open if args.fasta.endswith((".gz",".bgz")) else (bz2.open if args.fasta.endswith(".bz2") else open)

    with opener(args.fasta, "rt") as fa, open(args.out, "w") as gff:
        for chrom_raw, seq_one_hot in fasta_iter(fa):
            chrom = f"{args.id_prefix}{chrom_raw}"
            batch_offsets: List[int] = []
            batch_chunks : List[torch.Tensor] = []
            for off, chunk in sliding_windows(seq_one_hot):
                batch_offsets.append(off)
                batch_chunks.append(torch.from_numpy(chunk).float())
                if len(batch_chunks) == args.batch:
                    _run_batch(model, batch_chunks, batch_offsets, chrom, args.thr, gff)
                    batch_offsets, batch_chunks = [], []
            if batch_chunks:
                _run_batch(model, batch_chunks, batch_offsets, chrom, args.thr, gff)

    print(f"[DeepSplice] Wrote {args.out}")

if __name__ == "__main__":
    main()

# ───────────────────── Optional: checkpoint key converter ───────────────
# Save as fix_ckpt_keys.py to run once:
#   python fix_ckpt_keys.py src.pt dst.pt
if False:  # prevent execution when imported
    import argparse as _ap
    def _convert_main():
        ap = _ap.ArgumentParser("Rewrite DeepSplice checkpoint keys")
        ap.add_argument("src")
        ap.add_argument("dst")
        a = ap.parse_args()
        raw = _safe_torch_load(Path(a.src), device="cpu")
        state = _extract_state_dict(raw)
        dummy = OSAI()
        fixed = _rename_legacy_keys(state, dummy.state_dict().keys())
        good = {k:v for k,v in fixed.items() if k in dummy.state_dict() and dummy.state_dict()[k].shape == v.shape}
        torch.save(good, a.dst)
        print("Wrote", a.dst)
    _convert_main()
