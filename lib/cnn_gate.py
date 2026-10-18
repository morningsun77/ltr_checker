#!/usr/bin/env python
# coding: utf-8
"""
Standalone CNN detection stage of ltr_checker.

Runs ONLY the CNN sliding-window detector (model/model_cnn.pth) over a genome and
writes the predicted candidate LTR-RT regions (BED). No structural annotation, no
6-module filtering. GPU-accelerated with a vectorised one-hot encoder.

Output:
  <out>/<name>_cnn_regions.bed   chr  start  end
  <out>/<name>_cnn_summary.tsv   per-chromosome region count / bp / fraction
"""
import argparse, os, time
import numpy as np
import torch
import torch.nn.functional as F
from model import DETECT
from sequence_utils import genome_dict, ltr_range

# ASCII -> code (A,C,G,T,N -> 0..4 ; everything else -> 4/N), matches model.fasta2one_hot
LUT = np.full(256, 4, dtype=np.int64)
for ch, code in zip("ACGTN", range(5)):
    LUT[ord(ch)] = code
    LUT[ord(ch.lower())] = code

def encode_chrom(seq):
    return LUT[np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)]

def predict_chrom(codes, model, device, stride, win, batch):
    """Return a flat list of 0/1 window predictions along the chromosome."""
    n = len(codes)
    starts = list(range(0, n - win + 1, stride))  # only full windows, matches original (end<len)
    preds = []
    with torch.no_grad():
        for b0 in range(0, len(starts), batch):
            idx = starts[b0:b0 + batch]
            mat = np.stack([codes[s:s + win] for s in idx])         # (B, win)
            t = torch.from_numpy(mat).to(device)
            oh = F.one_hot(t, num_classes=5).permute(0, 2, 1).float()  # (B,5,win)
            oh = oh.unsqueeze(1)                                     # (B,1,5,win)
            out = model(oh)
            preds.extend(out.max(1)[1].cpu().numpy().tolist())
    return preds

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", required=True, help="output prefix, e.g. rice / maize")
    ap.add_argument("--model", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "model", "model_cnn.pth"))
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--stride", type=int, default=10000)
    ap.add_argument("--win", type=int, default=50000)
    ap.add_argument("--batch", type=int, default=128)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    device = args.device
    if device == "cuda" and not torch.cuda.is_available():
        print("[WARN] cuda requested but not available -> falling back to cpu", flush=True)
        device = "cpu"

    model = DETECT()
    model.load_state_dict(torch.load(args.model, map_location="cpu"))
    model.eval().to(device)

    if device == "cuda":  # verify the GPU arch is supported by this torch build; else fall back
        try:
            with torch.no_grad():
                model(torch.zeros(1, 1, 5, args.win, device=device))
            print(f"[INFO] GPU OK: {torch.cuda.get_device_name(0)}", flush=True)
        except Exception as e:
            print(f"[WARN] GPU forward failed ({e}); falling back to cpu", flush=True)
            device = "cpu"; model.to("cpu")

    genome = genome_dict(args.genome)
    print(f"[INFO] {args.name}: {len(genome)} sequences, device={device}", flush=True)

    bed = os.path.join(args.out, f"{args.name}_cnn_regions.bed")
    summ = os.path.join(args.out, f"{args.name}_cnn_summary.tsv")
    t0 = time.time()
    total_regions = 0
    with open(bed, "w") as fb, open(summ, "w") as fs:
        fs.write("chrom\tlength\tn_regions\tregion_bp\tfraction\ttime_s\n")
        for chrom, seq in genome.items():
            ct = time.time()
            L = len(seq)
            if L < args.win:
                fs.write(f"{chrom}\t{L}\t0\t0\t0.0000\t0.0\n"); continue
            codes = encode_chrom(seq)
            preds = predict_chrom(codes, model, device, args.stride, args.win, args.batch)
            regions = ltr_range(chrom, [preds], L, args.stride)  # list of "chr\tstart\tend"
            rbp = 0
            for r in regions:
                c, s, e = r.split("\t")
                s, e = int(s), min(int(e), L)
                if e > s:
                    fb.write(f"{c}\t{s}\t{e}\n"); rbp += e - s
            total_regions += len(regions)
            fb.flush()
            fs.write(f"{chrom}\t{L}\t{len(regions)}\t{rbp}\t{rbp/L:.4f}\t{time.time()-ct:.1f}\n"); fs.flush()
            print(f"[INFO] {chrom}: len={L} regions={len(regions)} bp={rbp} ({100*rbp/L:.1f}%) {time.time()-ct:.1f}s", flush=True)
    print(f"[DONE] {args.name}: {total_regions} candidate regions in {time.time()-t0:.1f}s -> {bed}", flush=True)

if __name__ == "__main__":
    main()
