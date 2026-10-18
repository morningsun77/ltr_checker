#!/usr/bin/env python3
"""
Normalize a tool's raw LTR-RT candidates to LTRharvest .scn format in GENOME coordinates,
so they can be fed to LTR_retriever -inharvest. Region/block-local coords are mapped to the
genome via the run's input.map.tsv (record id / index -> chrom, offset), and seq-nr is set to
the chromosome's index in a fixed chrom order.

LTRharvest .scn columns (1 line per element):
  s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim seq-nr
Coordinates are 1-based, per-chromosome (== genome position within that chromosome).
Only full elements (both LTRs present) are emitted.
"""
import argparse, os, sys, glob

def load_map(mapf):
    by_id, by_idx = {}, {}
    with open(mapf) as fh:
        next(fh)
        for line in fh:
            rid, idx, chrom, off, ln = line.rstrip("\n").split("\t")
            by_id[rid] = (chrom, int(off)); by_idx[int(idx)] = (chrom, int(off))
    return by_id, by_idx

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True, choices=["finder", "harvest", "detector", "look4ltrs"])
    ap.add_argument("--raw", required=True, help="raw file or run dir (detector/look4ltrs)")
    ap.add_argument("--map", required=True)
    ap.add_argument("--chrom-order", required=True, help="chroms.txt (order = seq-nr)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    by_id, by_idx = load_map(a.map)
    order = {c.strip(): i for i, c in enumerate(open(a.chrom_order)) if c.strip()}
    rows = []   # (chrom_idx, s, e, ls, le, rs, re, sim)

    def add(chrom, off, es, ee, ls, le, rs, re, sim):
        if chrom not in order:
            return
        ci = order[chrom]
        # to genome (within-chrom) 1-based
        es, ee = off + es, off + ee
        ls, le = off + ls, off + le
        rs, re = off + rs, off + re
        rows.append((ci, es, ee, ls, le, rs, re, sim))

    if a.tool == "finder":
        with open(a.raw) as fh:
            for line in fh:
                if not line.startswith("["):
                    continue
                c = line.rstrip("\n").split("\t")
                rid = c[1].split()[0]
                if rid not in by_id:
                    continue
                chrom, off = by_id[rid]
                s, e = [int(x) for x in c[2].split("-")]
                ll, rl = [int(x) for x in c[3].split(",")[:2]]
                sim = 100.0
                try: sim = float(c[-1]) * (100 if float(c[-1]) <= 1 else 1)
                except Exception: pass
                add(chrom, off, s, e, s, s + ll - 1, e - rl + 1, e, sim)

    elif a.tool == "harvest":
        with open(a.raw) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split()
                if len(f) < 11:
                    continue
                seqnr = int(f[10])
                if seqnr not in by_idx:
                    continue
                chrom, off = by_idx[seqnr]
                add(chrom, off, int(f[0]), int(f[1]), int(f[3]), int(f[4]),
                    int(f[6]), int(f[7]), float(f[9]))

    elif a.tool == "detector":
        beds = glob.glob(os.path.join(a.raw, "*Detector.bed")) if os.path.isdir(a.raw) else [a.raw]
        for bed in beds:
            for line in open(bed):
                c = line.rstrip("\n").split("\t")
                if not c or c[0] in ("SeqID", "") or c[0].startswith("Start") or len(c) < 7:
                    continue
                rid = c[0]
                if rid not in by_id:
                    continue
                chrom, off = by_id[rid]
                try:
                    es, ee, ls, le, rs, re = (int(c[1]), int(c[2]), int(c[3]),
                                              int(c[4]), int(c[5]), int(c[6]))
                    sim = float(c[7]) if c[7].replace('.', '').isdigit() else 90.0
                except ValueError:
                    continue
                add(chrom, off, es, ee, ls, le, rs, re, sim)

    elif a.tool == "look4ltrs":
        rtrs = glob.glob(os.path.join(a.raw, "Rtr", "*.rtr")) if os.path.isdir(a.raw) else [a.raw]
        for rtr in rtrs:
            for line in open(rtr):
                c = line.rstrip("\n").split("\t")
                if not c or c[0] in ("chrom", "") or len(c) < 6:
                    continue
                rid = c[0]
                if rid not in by_id or "NA" in (c[2], c[3], c[4], c[5]):
                    continue
                chrom, off = by_id[rid]
                try:
                    ls, le, rs, re = int(c[2]), int(c[3]), int(c[4]), int(c[5])
                    sim = float(c[15]) * 100 if len(c) > 15 and c[15] not in ("NA", "") else 90.0
                except ValueError:
                    continue
                add(chrom, off, ls, re, ls, le, rs, re, sim)   # element = LeftStart..RightEnd

    # emit sorted, valid (s<e, ltrs ordered)
    good = []
    for ci, es, ee, ls, le, rs, re, sim in rows:
        if es < ee and ls <= le < rs <= re:
            good.append((ci, es, ee, ls, le, rs, re, sim))
    good.sort()
    with open(a.out, "w") as oh:
        oh.write("# normalized to LTRharvest .scn (genome coords) from tool=%s\n" % a.tool)
        oh.write("# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim seq-nr\n")
        for ci, es, ee, ls, le, rs, re, sim in good:
            oh.write(f"{es} {ee} {ee-es+1} {ls} {le} {le-ls+1} {rs} {re} {re-rs+1} {sim:.2f} {ci}\n")
    print(f"[to_scn] {a.tool}: {len(good)} elements -> {a.out}")

if __name__ == "__main__":
    main()
