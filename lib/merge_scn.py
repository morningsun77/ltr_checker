#!/usr/bin/env python
# coding: utf-8
"""Overlap-aware, priority-resolved merge of per-tool .scn candidate files.

This replaces the old exact-coordinate dedup. Candidate calls from different tools that
hit the SAME locus (reciprocal overlap >= `reciprocal`, default 0.80) collapse to a
single representative: the call from the highest-priority tool. Tool priority is the
order the .scn files are given (first = most trustworthy); ltr_checker passes them as
LTR_Finder > LTRharvest > LtrDetector. Ties (same priority) break to the longer element,
then the smaller start, so the result is deterministic.

Recall safety: a call that has no reciprocal-overlap partner forms its own singleton
cluster and is always kept. The merge therefore only ever collapses near-duplicate
calls of one locus; it never drops a locus that only one tool found (e.g. LtrDetector's
unique additions) and so cannot reduce locus-level recall by construction. (The pipeline
additionally verifies this empirically against the rice/maize truth sets.)

The .scn columns are the LTRharvest -tabout form, already remapped to genome coords:
    s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim seq-nr
so field 0 = element start, field 1 = element end, field -1 = seq-nr (chrom index).
"""
import argparse
import sys
from collections import defaultdict


def load_calls(scn_paths):
    """Read every non-comment call from the given .scn files, tagging each with the
    priority (= index of its file in `scn_paths`, lower is better). Returns a list of
    (seqnr, start, end, prio, line) with start <= end."""
    calls = []
    for prio, path in enumerate(scn_paths):
        try:
            fh = open(path)
        except OSError:
            continue
        with fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.split()
                try:
                    start, end = int(f[0]), int(f[1])
                    seqnr = f[-1]
                except (IndexError, ValueError):
                    continue
                if end < start:
                    start, end = end, start
                if not line.endswith("\n"):
                    line += "\n"
                calls.append((seqnr, start, end, prio, line))
    return calls


def merge_calls(calls, reciprocal=0.80):
    """Collapse reciprocally-overlapping calls to one representative each.

    Returns (rep_lines, stats) where stats = {in, out, collapsed}.
    Greedy single-linkage against each cluster's current representative interval:
    processing calls sorted by start, a call joins the first active cluster it overlaps
    reciprocally, and promotes itself to representative if it is higher priority (or a
    better tie-break). Active clusters whose representative ends before the current call
    cannot overlap any later call (calls are start-sorted) and are retired from the scan.
    """
    by_seq = defaultdict(list)
    for c in calls:
        by_seq[c[0]].append(c)

    rep_lines = []
    n_in = len(calls)
    n_out = 0
    for seqnr, items in by_seq.items():
        items.sort(key=lambda c: (c[1], c[2]))
        clusters = []        # each: [start, end, prio, length, line]
        active = []          # indices into `clusters` that may still overlap
        for (_sq, s, e, prio, line) in items:
            length = e - s + 1
            # retire clusters that end before this call starts (sorted by start => they
            # cannot overlap this or any later call)
            active = [ci for ci in active if clusters[ci][1] >= s]
            placed = False
            for ci in active:
                cs, ce, cprio, clen, _cline = clusters[ci]
                ov = min(e, ce) - max(s, cs) + 1
                if ov <= 0:
                    continue
                if ov >= reciprocal * length and ov >= reciprocal * clen:
                    better = (prio < cprio) or (
                        prio == cprio and (length > clen or (length == clen and s < cs)))
                    if better:
                        clusters[ci] = [s, e, prio, length, line]
                    placed = True
                    break
            if not placed:
                clusters.append([s, e, prio, length, line])
                active.append(len(clusters) - 1)
        for cs, ce, cprio, clen, cline in clusters:
            rep_lines.append(cline)
        n_out += len(clusters)
    return rep_lines, {"in": n_in, "out": n_out, "collapsed": n_in - n_out}


def write_merged(scn_paths, out_path, reciprocal=0.80):
    """Merge `scn_paths` (priority order) and write the representatives to `out_path`.
    Returns the stats dict."""
    calls = load_calls(scn_paths)
    rep_lines, stats = merge_calls(calls, reciprocal=reciprocal)
    with open(out_path, "w") as oh:
        oh.writelines(rep_lines)
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scn", nargs="+", required=True,
                    help="per-tool .scn files in PRIORITY order (first = most trustworthy)")
    ap.add_argument("--out", required=True, help="merged candidates .scn")
    ap.add_argument("--reciprocal", type=float, default=0.80,
                    help="reciprocal-overlap fraction to treat two calls as the same locus")
    a = ap.parse_args()
    stats = write_merged(a.scn, a.out, reciprocal=a.reciprocal)
    print(f"[merge_scn] in={stats['in']} out={stats['out']} "
          f"collapsed={stats['collapsed']} -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
