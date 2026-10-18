#!/usr/bin/env python3
"""
Parallel candidate-finder runners for the ltr_checker `candidates` stage.

Input  : <candidates>/input_dir/input.fa  (region pieces named s{idx}, built by
         make_input.py) + input.map.tsv.
Each finder saturates its OWN thread budget. The driver launches the three finders
either one-at-a-time (serial, default) or concurrently (--parallel-tools).

Every runner leaves its raw output in the exact shape lib/to_scn.py already expects:
  finder   -> <td>/finder.w2.scn  (ltr_finder -w 2 '[' table rows, keyed by s{idx};
              to_scn reads these rows directly)
  harvest  -> <td>/harvest.scn  (LTRharvest .scn whose seq-nr column has been rewritten
              from each worker's LOCAL index to the GLOBAL record index)
  detector -> <td>  (LtrDetector destDir; to_scn globs *Detector.bed there)

How each tool reaches its thread budget (all three use the same bundle strategy — pieces
greedily balanced into N bundles, one worker process per bundle):
  ltr_finder   raw ltr_finder (-w 2) run once per bundle. We drive the ltr_finder binary
               directly (NOT the LTR_FINDER_parallel perl wrapper) so the CLI parameters
               apply; the piece-splitting + parallel strategy mirrors LTR_FINDER_parallel.
  ltr_harvest  one gt suffixerator index + gt ltrharvest per bundle (both single-threaded).
  ltrdetector  LtrDetector's own -nThreads (single run).
"""
import os, sys, glob, subprocess, multiprocessing as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from faio import read_fasta, write_fasta


# --------------------------------------------------------------------------- helpers
def _read_pieces(indir):
    """input.fa -> [(record_id, seq), ...] in write order."""
    fa = os.path.join(indir, "input.fa")
    return [(name.split()[0], seq) for name, seq in read_fasta(fa)]


def _bundle(records, n):
    """Greedy longest-processing-time packing of records into <= n balanced bundles."""
    if not records:
        return []
    n = max(1, min(n, len(records)))
    load = [0] * n
    bins = [[] for _ in range(n)]
    for rid, seq in sorted(records, key=lambda r: len(r[1]), reverse=True):
        j = min(range(n), key=lambda k: load[k])
        bins[j].append((rid, seq))
        load[j] += len(seq)
    return [b for b in bins if b]


# --------------------------------------------------------------------------- ltr_finder
def _finder_worker(arg):
    ltr_finder, bundle, td, p, bi = arg
    fa = os.path.join(td, f"b{bi}.fa")
    with open(fa, "w") as fh:
        for gid, seq in bundle:
            write_fasta(fh, gid, seq)
    # -w 2 table output; -M is a fraction (identity/100). Records keep their s{idx} name,
    # which to_scn maps back to genome coordinates, so no coordinate stitching is needed.
    cmd = (f'"{ltr_finder}" -w 2 -C -D {p["max"]} -d {p["min"]} '
           f'-L {p["max_ltr"]} -l {p["min_ltr"]} -p 20 -M {p["identity"] / 100.0:.2f} "{fa}"')
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=td)
    with open(os.path.join(td, f"b{bi}.log"), "w") as lh:
        lh.write(out.stderr or "")
    return [ln for ln in out.stdout.splitlines() if ln.startswith("[")]


def run_finder(ltr_finder, indir, td, p, nthreads):
    """Parallel raw ltr_finder: bundle the region pieces and run one ltr_finder -w 2 per
    bundle, collecting the '[' table rows (which to_scn reads directly)."""
    os.makedirs(td, exist_ok=True)
    bundles = _bundle(_read_pieces(indir), nthreads)
    raw = os.path.join(td, "finder.w2.scn")
    if not bundles:
        open(raw, "w").close()
        return raw
    args = [(ltr_finder, b, td, p, i) for i, b in enumerate(bundles)]
    if len(args) == 1:
        results = [_finder_worker(args[0])]
    else:
        with mp.Pool(processes=len(args)) as pool:
            results = pool.map(_finder_worker, args)
    with open(raw, "w") as oh:
        for r in results:
            oh.write("\n".join(r) + ("\n" if r else ""))
    return raw


# --------------------------------------------------------------------------- ltr_harvest
def _harvest_worker(arg):
    gt, bundle, td, p, bi = arg
    fa = os.path.join(td, f"b{bi}.fa")
    gids = []
    with open(fa, "w") as fh:
        for gid, seq in bundle:
            write_fasta(fh, gid, seq)
            gids.append(gid)                      # local position -> global record id
    idxname = os.path.join(td, f"b{bi}.idx")
    scn = os.path.join(td, f"b{bi}.scn")
    log = os.path.join(td, f"b{bi}.log")
    with open(log, "w") as lh:
        subprocess.run(
            f'"{gt}" suffixerator -db "{fa}" -indexname "{idxname}" '
            f'-tis -suf -lcp -des -ssp -sds -dna',
            shell=True, check=True, stdout=lh, stderr=lh)
        subprocess.run(
            f'"{gt}" ltrharvest -index "{idxname}" '
            f'-minlenltr {p["min_ltr"]} -maxlenltr {p["max_ltr"]} '
            f'-mindistltr {p["min"]} -maxdistltr {p["max"]} '
            f'-similar {p["identity"]} -vic 10 -seed 20 -seqids yes > "{scn}"',
            shell=True, check=True, stderr=lh)
    # rewrite column 11 (local seq-nr) -> global record index, so to_scn's by_idx works
    out = []
    for line in open(scn):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 11:
            continue
        local = int(f[10])
        if local >= len(gids):
            continue
        f[10] = gids[local][1:]                   # 's{G}' -> 'G'
        out.append(" ".join(f))
    return out


def run_harvest(gt, indir, td, p, nthreads):
    """Parallel LTRharvest: bundle the region pieces and run one gt per bundle."""
    os.makedirs(td, exist_ok=True)
    bundles = _bundle(_read_pieces(indir), nthreads)
    if not bundles:
        raw = os.path.join(td, "harvest.scn")
        open(raw, "w").close()
        return raw
    args = [(gt, b, td, p, i) for i, b in enumerate(bundles)]
    if len(args) == 1:
        results = [_harvest_worker(args[0])]
    else:
        with mp.Pool(processes=len(args)) as pool:
            results = pool.map(_harvest_worker, args)
    raw = os.path.join(td, "harvest.scn")
    with open(raw, "w") as oh:
        for r in results:
            oh.write("\n".join(r) + ("\n" if r else ""))
    return raw


# --------------------------------------------------------------------------- ltrdetector
def run_detector(ltrdetector, indir, td, p, nthreads):
    """Run LtrDetector on a private copy of input.fa (avoids racing other tools' files)."""
    os.makedirs(td, exist_ok=True)
    det_in = os.path.join(td, "in")
    os.makedirs(det_in, exist_ok=True)
    link = os.path.join(det_in, "input.fa")
    if not os.path.exists(link):
        os.symlink(os.path.abspath(os.path.join(indir, "input.fa")), link)
    log = os.path.join(td, "tool.log")
    cmd = (f'"{ltrdetector}" -fasta "{det_in}" -destDir "{td}" '
           f'-minLen {p["min"]} -maxLen {p["max"]} '
           f'-minLenLTR {p["min_ltr"]} -maxLenLTR {p["max_ltr"]} '
           f'-id {p["identity"]} -nThreads {nthreads} > "{log}" 2>&1')
    subprocess.run(cmd, shell=True, check=True)
    return td
