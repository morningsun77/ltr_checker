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
  ltr_harvest  one gt suffixerator index + gt ltrharvest per bundle (both single-threaded);
               each gt run is capped by HARVEST_TIMEOUT and, on timeout, the bundle is
               retried piece-by-piece so one stalling piece can't hang the whole run.
  ltrdetector  LtrDetector's own -nThreads (single run).
"""
import os, sys, glob, signal, subprocess, multiprocessing as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from faio import read_fasta, write_fasta

# LTRharvest can grow worse-than-linearly on repeat-dense sequence and occasionally stalls
# for hours on a single pathological piece. We cap each gt run with a timeout (seconds) and,
# on timeout, retry the bundle piece-by-piece, skipping only the piece that still stalls.
# 0 disables the cap. Override with the LTRC_HARVEST_TIMEOUT environment variable.
HARVEST_TIMEOUT = int(os.environ.get("LTRC_HARVEST_TIMEOUT", "1800"))


class _Timeout(Exception):
    pass


def _run_timed(cmd, log, timeout):
    """Run a shell command, appending output to `log`. Kill the whole process group on
    timeout (so the gt child dies too) and raise _Timeout. timeout<=0 means no cap."""
    with open(log, "a") as lh:
        proc = subprocess.Popen(cmd, shell=True, stdout=lh, stderr=lh,
                                start_new_session=True)
        try:
            rc = proc.wait(timeout=timeout if timeout and timeout > 0 else None)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass
            proc.wait()
            raise _Timeout(cmd)
        if rc != 0:
            raise subprocess.CalledProcessError(rc, cmd)


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
    # argv list (no shell): every path/parameter is a separate argument, so a space or
    # shell metacharacter in a path can neither break the command nor inject into it.
    cmd = [ltr_finder, "-w", "2", "-C", "-D", str(p["max"]), "-d", str(p["min"]),
           "-L", str(p["max_ltr"]), "-l", str(p["min_ltr"]),
           "-p", str(p.get("finder_p", 20)),
           "-M", f"{p['identity'] / 100.0:.2f}"]
    trna = p.get("trna")
    if trna:                       # optional PBS tRNA database (ltr_finder -s FASTA)
        cmd += ["-s", trna]
    cmd += [fa]
    out = subprocess.run(cmd, capture_output=True, text=True, cwd=td)
    with open(os.path.join(td, f"b{bi}.log"), "w") as lh:
        lh.write(out.stderr or "")
    # Fail loud: a nonzero return (e.g. a missing/incompatible binary -> rc 127) must NOT
    # be mistaken for "no LTR-RTs found" and silently yield an empty candidate set.
    if out.returncode != 0:
        raise RuntimeError(
            f"ltr_finder failed (rc={out.returncode}) on bundle b{bi}; see "
            f"{os.path.join(td, f'b{bi}.log')}. stderr tail: "
            f"{(out.stderr or '').strip()[-300:]}")
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
def _harvest_group(gt, group, td, p, tag, timeout):
    """Run gt suffixerator + ltrharvest over `group` (a list of (gid, seq)); return .scn
    lines with column 11 rewritten from the local seq-nr to the global record id. Raise
    _Timeout if the (capped) ltrharvest run exceeds `timeout`."""
    fa = os.path.join(td, f"{tag}.fa")
    gids = []
    with open(fa, "w") as fh:
        for gid, seq in group:
            write_fasta(fh, gid, seq)
            gids.append(gid)                      # local position -> global record id
    idxname = os.path.join(td, f"{tag}.idx")
    scn = os.path.join(td, f"{tag}.scn")
    log = os.path.join(td, f"{tag}.log")
    open(log, "w").close()
    # suffixerator (indexing) is bounded; the cap is really for ltrharvest below.
    _run_timed(f'"{gt}" suffixerator -db "{fa}" -indexname "{idxname}" '
               f'-tis -suf -lcp -des -ssp -sds -dna', log, timeout)
    _run_timed(f'"{gt}" ltrharvest -index "{idxname}" '
               f'-minlenltr {p["min_ltr"]} -maxlenltr {p["max_ltr"]} '
               f'-mindistltr {p["min"]} -maxdistltr {p["max"]} '
               f'-similar {p["identity"]} -vic {p.get("vic", 10)} -seed {p.get("seed", 20)} '
               f'-seqids yes > "{scn}"',
               log, timeout)
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


def _harvest_worker(arg):
    gt, bundle, td, p, bi = arg
    try:                                          # fast path: whole bundle in one index
        return _harvest_group(gt, bundle, td, p, f"b{bi}", HARVEST_TIMEOUT)
    except _Timeout:
        # a piece in this bundle stalled ltrharvest -> retry piece-by-piece and skip only
        # the one that still exceeds the cap, so one bad piece can't hang the whole run.
        sys.stderr.write(f"[run_finders] harvest bundle b{bi} timed out "
                         f"(>{HARVEST_TIMEOUT}s); retrying its {len(bundle)} pieces singly\n")
        out = []
        for j, rec in enumerate(bundle):
            try:
                out.extend(_harvest_group(gt, [rec], td, p, f"b{bi}_p{j}", HARVEST_TIMEOUT))
            except _Timeout:
                sys.stderr.write(f"[run_finders] harvest SKIPPED piece {rec[0]} "
                                 f"(ltrharvest >{HARVEST_TIMEOUT}s on a single piece)\n")
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
    # argv list (no shell); stdout+stderr captured to the log file. check=True already
    # fails loud on a nonzero return, so a broken LtrDetector cannot pass silently.
    cmd = [ltrdetector, "-fasta", det_in, "-destDir", td,
           "-minLen", str(p["min"]), "-maxLen", str(p["max"]),
           "-minLenLTR", str(p["min_ltr"]), "-maxLenLTR", str(p["max_ltr"]),
           "-id", str(p["identity"]), "-nThreads", str(nthreads)]
    with open(log, "w") as lh:
        subprocess.run(cmd, stdout=lh, stderr=subprocess.STDOUT, check=True)
    return td
