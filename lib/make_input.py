#!/usr/bin/env python3
"""
Build the tool-input FASTA + coordinate map for one benchmark run.

Two modes:
  alone : one record per (selected) chromosome  -> record covers whole chrom (offset 0)
  gate  : one record per CNN candidate region   -> record = chrom[start:end] (offset=start)

Every record is named s{index} (index = write order) so that all four downstream tools
(LFP/ltrharvest/LtrDetector/Look4LTRs) can be remapped to genome coordinates uniformly:
  - LFP / LtrDetector / Look4LTRs report the record NAME  -> look up by id 's{index}'
  - ltrharvest reports the record INDEX (seq-nr)          -> look up by index

Any record longer than --chunk bp is cut into overlapping pieces of at most --chunk bp
(overlap --overlap bp), each written as its own s{index} record with offset advanced by
the piece start. This lets the finders run on bounded pieces in parallel; the overlap
must exceed the largest possible element so none is lost on a cut line (the driver sizes
it from the LTR-separation / LTR-length limits). Overlap-zone duplicate calls are deduped
downstream when the per-tool .scn files are merged.

Outputs:
  <outdir>/input_dir/input.fa      multi-FASTA (dir form needed by LtrDetector/Look4LTRs)
  <outdir>/input.map.tsv           id  index  chrom  offset  length
"""
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from faio import read_fasta, write_fasta


def pieces(seq, chunk, overlap):
    """Yield (piece_start, piece_seq): whole seq if <= chunk, else overlapping windows."""
    n = len(seq)
    if chunk <= 0 or n <= chunk:
        yield 0, seq
        return
    step = max(1, chunk - overlap)
    start = 0
    while start < n:
        end = min(start + chunk, n)
        yield start, seq[start:end]
        if end >= n:
            break
        start += step

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True, help="genome (or single-chrom) FASTA")
    ap.add_argument("--mode", required=True, choices=["alone", "gate"])
    ap.add_argument("--regions", help="CNN regions BED (required for gate): chrom start end")
    ap.add_argument("--chroms", help="comma-separated chrom whitelist (default: all in genome)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-region", type=int, default=0,
                    help="drop gate regions shorter than this many bp (default 0 = keep all)")
    ap.add_argument("--chunk", type=int, default=5_000_000,
                    help="cut records longer than this into overlapping pieces (bp, 0=off)")
    ap.add_argument("--overlap", type=int, default=100_000,
                    help="overlap between pieces (bp); must exceed the largest element")
    args = ap.parse_args()

    keep = set(args.chroms.split(",")) if args.chroms else None
    indir = os.path.join(args.outdir, "input_dir")
    os.makedirs(indir, exist_ok=True)
    fa = os.path.join(indir, "input.fa")
    mapf = os.path.join(args.outdir, "input.map.tsv")

    # load only the chromosomes we need
    genome = {}
    for name, seq in read_fasta(args.genome):
        if keep is None or name in keep:
            genome[name] = seq
    if not genome:
        sys.exit(f"[make_input] ERROR: no chromosomes matched (keep={keep})")

    idx = 0
    with open(fa, "w") as fh, open(mapf, "w") as mh:
        mh.write("id\tindex\tchrom\toffset\tlength\n")

        def emit(chrom, base, seq):
            """Write seq (cut into pieces) as s{index} records; offset tracks genome pos."""
            nonlocal idx
            for p_start, p_seq in pieces(seq, args.chunk, args.overlap):
                rid = f"s{idx}"
                write_fasta(fh, rid, p_seq)
                mh.write(f"{rid}\t{idx}\t{chrom}\t{base + p_start}\t{len(p_seq)}\n")
                idx += 1

        if args.mode == "alone":
            for name in genome:
                emit(name, 0, genome[name])
        else:  # gate
            if not args.regions:
                sys.exit("[make_input] ERROR: --regions required for gate mode")
            n_skip_chr = n_skip_len = n_reg = 0
            with open(args.regions) as rf:
                for line in rf:
                    if not line.strip() or line.startswith(("#", "track", "chrom\t")):
                        continue
                    f = line.split()
                    chrom, s, e = f[0], int(f[1]), int(f[2])
                    if chrom not in genome:
                        n_skip_chr += 1; continue
                    s = max(0, s); e = min(e, len(genome[chrom]))
                    if e - s < args.min_region or e <= s:
                        n_skip_len += 1; continue
                    emit(chrom, s, genome[chrom][s:e])
                    n_reg += 1
            sys.stderr.write(f"[make_input] gate: {n_reg} regions -> {idx} pieces, "
                             f"skipped {n_skip_chr} off-target chrom, {n_skip_len} short/empty\n")
    print(f"[make_input] mode={args.mode} records={idx} -> {fa}")

if __name__ == "__main__":
    main()
