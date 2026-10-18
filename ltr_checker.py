#!/usr/bin/env python
# coding: utf-8
"""
ltr_checker - a CNN-gated structural pipeline for LTR retrotransposon annotation.

Pipeline (four stages):
  gate       A CNN scans the genome in 50 kb sliding windows and flags the regions
             likely to contain LTR-RTs  ->  <out>/gate/<name>_cnn_regions.bed
  candidates Three structural tools (LTR_Finder, LTRharvest, LtrDetector)
             run ONLY inside the flagged regions; their calls are remapped to genome
             coordinates and merged  ->  <out>/candidates/candidates.scn
  intact     LTR_retriever adjudicates the merged candidates into intact elements and
             builds a RepeatMasker library  ->  <out>/intact/genome.fa.pass.list
                                                 <out>/intact/genome.fa.LTRlib.fa
  annotate   RepeatMasker annotates the whole genome with that library, using either
             the EDTA or the HiTE parameter preset  ->  <out>/annotate/

ltr_checker's core job is the CANDIDATE SET (gate + candidates); that is what runs by
default. The intact (LTR_retriever) and annotate (RepeatMasker) stages are optional
downstream steps, enabled with --stage intact,annotate (or --stage all).

External tools are located through config.yaml (or the same-named environment
variables); see README.md.

Author: Zhao-Yang Chen
"""
import argparse, os, sys, subprocess, shutil, time, glob

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "lib")
STAGES = ["gate", "candidates", "intact", "annotate"]
# to_scn tool names differ from the CLI/binary names
SCN_TOOL = {"ltr_finder": "finder", "ltr_harvest": "harvest", "ltrdetector": "detector"}


def log(msg):
    print(f"[ltr_checker] {msg}", flush=True)


def sh(cmd, **kw):
    """Run a shell command, stream output, raise on failure."""
    log("$ " + (cmd if isinstance(cmd, str) else " ".join(cmd)))
    subprocess.run(cmd, shell=isinstance(cmd, str), check=True, **kw)


def load_config(path):
    """Minimal `key: value` reader (no PyYAML dependency). Environment variables of
    the same name win; missing keys fall back to the bare binary name (found via PATH)."""
    defaults = {
        # ltr_finder (the raw LTR_Finder binary, bundled inside the LTR_FINDER_parallel
        # package) and LtrDetector ship bundled under tools/ (resolved against the package
        # root below). gt / LTR_retriever / RepeatMasker come from the conda environment
        # and are looked up on $PATH.
        "ltr_finder": "tools/LTR_FINDER_parallel/bin/LTR_FINDER.x86_64-1.0.7/ltr_finder",
        "gt": "gt",
        "LtrDetector": "tools/LtrDetector/bin/LtrDetector",
        "LTR_retriever": "LTR_retriever",
        "RepeatMasker": "RepeatMasker",
    }
    cfg = dict(defaults)
    if path and os.path.exists(path):
        for line in open(path):
            line = line.split("#", 1)[0].strip()
            if not line or ":" not in line:
                continue
            k, v = line.split(":", 1)
            k, v = k.strip(), v.strip().strip('"').strip("'")
            if v:
                cfg[k] = v
    for k in cfg:
        if os.environ.get(k):
            cfg[k] = os.environ[k]
    # Resolve bundled-tool paths (relative and containing a path separator) against
    # the package root, so they work regardless of the current working directory
    # (several stages cd into other dirs before invoking a tool). Bare binary names
    # with no separator are left alone and found via $PATH.
    for k, v in cfg.items():
        if v and not os.path.isabs(v) and os.sep in v:
            cand = os.path.join(HERE, v)
            if os.path.exists(cand):
                cfg[k] = os.path.abspath(cand)
    return cfg


def genome_chrom_order(genome, out):
    """Write chroms.txt = FASTA record order (used as the .scn seq-nr index)."""
    order = os.path.join(out, "chroms.txt")
    with open(order, "w") as oh:
        for line in open(genome):
            if line.startswith(">"):
                oh.write(line[1:].split()[0] + "\n")
    return order


# --------------------------------------------------------------------------- stages
def stage_gate(a, cfg):
    d = os.path.join(a.output, "gate")
    os.makedirs(d, exist_ok=True)
    sh([sys.executable, os.path.join(LIB, "cnn_gate.py"),
        "--genome", a.genome, "--out", d, "--name", a.name,
        "--model", a.model, "--device", a.device,
        "--stride", str(a.stride), "--win", "50000"], cwd=LIB)
    bed = os.path.join(d, f"{a.name}_cnn_regions.bed")
    if not os.path.exists(bed):
        sys.exit("gate stage produced no regions BED")
    log(f"gate -> {bed}")
    return bed


def stage_candidates(a, cfg, regions):
    import concurrent.futures as cf
    sys.path.insert(0, LIB)
    import run_finders as rf
    d = os.path.join(a.output, "candidates")
    os.makedirs(d, exist_ok=True)
    indir = os.path.join(d, "input_dir")
    mapf = os.path.join(d, "input.map.tsv")
    # overlap between 5 Mb pieces must exceed the largest possible element so none is
    # lost on a cut line: max element span = max LTR separation + the two LTR lengths.
    chunk = 5_000_000
    overlap = max(100_000, a.max + 2 * a.max_ltr)
    # 1) build the region sub-sequences (cut into overlapping pieces) shared by all tools
    sh([sys.executable, os.path.join(LIB, "make_input.py"),
        "--genome", a.genome, "--mode", "gate", "--regions", regions, "--outdir", d,
        "--chunk", str(chunk), "--overlap", str(overlap)], cwd=LIB)
    order = genome_chrom_order(a.genome, d)
    params = {"min": a.min, "max": a.max, "min_ltr": a.min_ltr,
              "max_ltr": a.max_ltr, "identity": a.identity}

    # 2) one callable per requested tool: run the finder, then normalise to genome .scn
    def run_tool(tool):
        td = os.path.join(d, tool)
        nt = {"ltr_finder": a.finder_threads, "ltr_harvest": a.harvest_threads,
              "ltrdetector": a.detector_threads}[tool]
        log(f"[{tool}] starting ({nt} threads)")
        if tool == "ltr_finder":
            raw = rf.run_finder(cfg["ltr_finder"], indir, td, params, nt)
        elif tool == "ltr_harvest":
            raw = rf.run_harvest(cfg["gt"], indir, td, params, nt)
        elif tool == "ltrdetector":
            raw = rf.run_detector(cfg["LtrDetector"], indir, td, params, nt)
        else:
            sys.exit(f"unknown tool {tool}")
        scn = os.path.join(d, f"{tool}.scn")
        sh([sys.executable, os.path.join(LIB, "to_scn.py"),
            "--tool", SCN_TOOL[tool], "--raw", raw, "--map", mapf,
            "--chrom-order", order, "--out", scn], cwd=LIB)
        log(f"[{tool}] done -> {scn}")
        return scn

    # 3) run the finders serially (default) or concurrently (--parallel-tools)
    if a.parallel_tools and len(a.tools) > 1:
        log(f"running {len(a.tools)} finders concurrently")
        with cf.ThreadPoolExecutor(max_workers=len(a.tools)) as ex:
            scns = list(ex.map(run_tool, a.tools))
    else:
        scns = [run_tool(t) for t in a.tools]

    # 4) merge the per-tool .scn files, deduping overlap-zone repeats
    merged = os.path.join(d, "candidates.scn")
    seen = set()
    n = 0
    with open(merged, "w") as oh:
        for s in scns:
            if not os.path.exists(s):
                continue
            for line in open(s):
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split()
                key = (f[-1], f[0], f[1])          # seq-nr, element start, element end
                if key in seen:
                    continue
                seen.add(key)
                oh.write(line)
                n += 1
    log(f"candidates -> {merged}  ({n} candidate LTR-RTs from {len(a.tools)} tools)")
    return merged


def _retriever_output(d, suffix):
    """Resolve an LTR_retriever output file in dir `d`.

    LTR_retriever standardises sequence names into a `genome.fa.mod` copy and then
    writes its outputs against that, so the real files carry a `.mod` infix
    (e.g. genome.fa.mod.LTRlib.fa). Match both the `.mod` and the plain form, and
    prefer the one that exists so the annotate stage finds the library it built.
    """
    for name in ("genome.fa.mod" + suffix, "genome.fa" + suffix):
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    # fall back to a glob, but skip the non-TGCA (.nmtf.) subset files
    hits = [h for h in sorted(glob.glob(os.path.join(d, "genome.fa*" + suffix)))
            if ".nmtf." not in os.path.basename(h)]
    return hits[0] if hits else os.path.join(d, "genome.fa" + suffix)


def stage_intact(a, cfg, scn):
    d = os.path.join(a.output, "intact")
    os.makedirs(d, exist_ok=True)
    g = os.path.join(d, "genome.fa")
    if not os.path.exists(g):
        os.symlink(os.path.abspath(a.genome), g)
    sh(f'cd "{d}" && "{cfg["LTR_retriever"]}" -genome genome.fa -inharvest "{scn}" '
       f'-threads {a.threads} 2>&1 | tail -60')
    passlist = _retriever_output(d, ".pass.list")
    lib_fa = _retriever_output(d, ".LTRlib.fa")
    if os.path.exists(passlist):
        nintact = sum(1 for l in open(passlist) if not l.startswith("#"))
        log(f"intact -> {passlist}  ({nintact} intact LTR-RTs)")
    else:
        log("intact -> no intact LTR-RTs passed (no pass.list); "
            "the candidate region may simply contain none")
    if os.path.exists(lib_fa):
        log(f"library -> {lib_fa}")
    return lib_fa


def stage_annotate(a, cfg, lib_fa):
    d = os.path.join(a.output, "annotate")
    os.makedirs(d, exist_ok=True)
    ref = os.path.join(d, "genome.fa")
    if not os.path.exists(ref):
        shutil.copy(a.genome, ref)
    pa = max(1, a.threads // 4)  # RepeatMasker -pa spawns 4 workers per unit
    if a.rm_preset == "edta":
        sh(f'"{cfg["RepeatMasker"]}" -e ncbi -pa {pa} -q -no_is -nolow -div 40 '
           f'-lib "{lib_fa}" "{ref}" 2>&1 | tail -30')
    else:  # hite
        sh(f'"{cfg["RepeatMasker"]}" -e ncbi -pa {pa} -q -no_is -nolow -div 40 -cutoff 225 -gff '
           f'-lib "{lib_fa}" "{ref}" 2>&1 | tail -30')
    out = ref + ".out"
    if os.path.exists(out):
        log(f"annotate ({a.rm_preset} preset) -> {out}")
    else:
        log("WARNING: RepeatMasker produced no .out")


# --------------------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(
        description="ltr_checker - CNN + 3 structural tools -> LTR-RT candidate set. "
                    "LTR_retriever (intact) and RepeatMasker (annotate) are optional "
                    "downstream stages, enabled with --stage.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # --- core: CNN gate + three tools -> candidate set -----------------------------
    core = p.add_argument_group(
        "candidate-set options (CNN gate + LTR_Finder / LTRharvest / LtrDetector)")
    core.add_argument("--genome", required=True, help="genome FASTA")
    core.add_argument("--output", required=True, help="output directory")
    core.add_argument("--stage", default="gate,candidates",
                      help="comma list of {gate,candidates,intact,annotate} or 'all'; "
                           "default runs the candidate set only (gate,candidates)")
    core.add_argument("--tools", default="ltr_finder,ltr_harvest,ltrdetector",
                      help="candidate tools to run (comma list), default all three")
    core.add_argument("--threads", type=int, default=10,
                      help="unified thread count for all tools; each finder uses it "
                           "unless overridden below")
    core.add_argument("--parallel-tools", action="store_true",
                      help="run the three candidate tools concurrently (peak threads = "
                           "sum of their counts) instead of one after another (default)")
    core.add_argument("--finder-threads", type=int, default=None,
                      help="parallel LTR_Finder workers (default: --threads)")
    core.add_argument("--harvest-threads", type=int, default=None,
                      help="parallel LTRharvest workers (default: --threads)")
    core.add_argument("--detector-threads", type=int, default=None,
                      help="threads for LtrDetector -nThreads (default: --threads)")
    core.add_argument("--device", default="cpu", choices=["cpu", "cuda"],
                      help="device for the CNN gate")
    core.add_argument("--model", default=os.path.join(HERE, "model", "model_cnn.pth"),
                      help="CNN weights")
    core.add_argument("--name", default="genome", help="output prefix")
    core.add_argument("--stride", type=int, default=10000, help="CNN window stride (bp)")
    # structural parameters shared by all three finders
    core.add_argument("--min", type=int, default=1000, help="min LTR separation (element internal), bp")
    core.add_argument("--max", type=int, default=15000, help="max LTR separation (element internal), bp")
    core.add_argument("--min-ltr", type=int, default=100, help="min LTR length, bp")
    core.add_argument("--max-ltr", type=int, default=7000, help="max LTR length, bp")
    core.add_argument("--identity", type=int, default=85, help="min LTR-pair identity (percent)")

    # --- optional downstream stages ------------------------------------------------
    down = p.add_argument_group(
        "optional downstream stages (enable with --stage intact,annotate or --stage all)")
    down.add_argument("--rm-preset", default="edta", choices=["edta", "hite"],
                      help="RepeatMasker preset for the annotate stage")

    # --- misc ----------------------------------------------------------------------
    p.add_argument("--config", default=os.path.join(HERE, "config.yaml"),
                   help="external-tool paths")
    a = p.parse_args()

    a.genome = os.path.abspath(a.genome)
    a.output = os.path.abspath(a.output)
    a.model = os.path.abspath(a.model)
    os.makedirs(a.output, exist_ok=True)
    want = STAGES if a.stage == "all" else [s.strip() for s in a.stage.split(",")]
    bad = [s for s in want if s not in STAGES]
    if bad:
        sys.exit(f"unknown stage(s): {bad}; choose from {STAGES}")
    a.tools = [t.strip() for t in a.tools.split(",")]
    # per-tool thread counts inherit --threads unless explicitly overridden
    a.finder_threads = a.finder_threads or a.threads
    a.harvest_threads = a.harvest_threads or a.threads
    a.detector_threads = a.detector_threads or a.threads
    cfg = load_config(a.config)

    log("=" * 60)
    log(f"genome={a.genome}")
    mode = "parallel" if a.parallel_tools else "serial"
    log(f"stages={want}  tools={a.tools}  rm_preset={a.rm_preset}")
    log(f"threads={a.threads}  finders={mode} "
        f"(finder={a.finder_threads} harvest={a.harvest_threads} detector={a.detector_threads})")
    log("=" * 60)
    t0 = time.time()

    regions = os.path.join(a.output, "gate", f"{a.name}_cnn_regions.bed")
    scn = os.path.join(a.output, "candidates", "candidates.scn")
    lib_fa = _retriever_output(os.path.join(a.output, "intact"), ".LTRlib.fa")

    if "gate" in want:
        regions = stage_gate(a, cfg)
    if "candidates" in want:
        if not os.path.exists(regions):
            sys.exit("candidates stage needs the gate BED; run --stage gate first")
        scn = stage_candidates(a, cfg, regions)
    if "intact" in want:
        if not os.path.exists(scn):
            sys.exit("intact stage needs candidates.scn; run --stage candidates first")
        lib_fa = stage_intact(a, cfg, scn)
    if "annotate" in want:
        if not os.path.exists(lib_fa):
            sys.exit("annotate stage needs the LTR library; run --stage intact first")
        stage_annotate(a, cfg, lib_fa)

    log(f"done in {time.time() - t0:.1f}s -> {a.output}")


if __name__ == "__main__":
    main()
