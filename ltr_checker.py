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
import argparse, os, sys, subprocess, shutil, time, glob, shlex

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
        # External tools. ltr_finder is on bioconda (v1.07) and gt / LTR_retriever /
        # RepeatMasker come from the conda environment, so all are bare names found on
        # $PATH. LtrDetector is NOT on bioconda; it ships bundled under tools/ and is
        # resolved to the bundled binary below when it is not on $PATH. ltr_finder has
        # the same bundled fallback for machines without the conda package.
        "ltr_finder": "ltr_finder",
        "gt": "gt",
        "LtrDetector": "LtrDetector",
        "LTR_retriever": "LTR_retriever",
        "RepeatMasker": "RepeatMasker",
        # Scientific parameters exposed here (rather than buried in the code) so they are
        # tunable and traceable. Values are strings in the file; the code casts as needed.
        "ltrharvest_vic": "10",      # LTRharvest -vic  (max TSD search radius, bp)
        "ltrharvest_seed": "20",     # LTRharvest -seed (min seed length for the LTR pair)
        "ltr_finder_p": "20",        # ltr_finder -p    (PBS/PPT search range)
        "ltr_finder_trna": "",       # optional ltr_finder -s tRNA FASTA (empty = none)
        "rm_div": "40",              # RepeatMasker -div    (max divergence, both presets)
        "rm_cutoff": "225",          # RepeatMasker -cutoff (hite preset only)
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
    # Resolve each tool path: an absolute path is used as-is; a relative path (containing a
    # separator, e.g. an explicit bundled path) is resolved against the package root so it
    # works regardless of the current working directory (several stages cd elsewhere before
    # invoking a tool); a BARE binary name is taken from $PATH, and for the tools we bundle
    # (ltr_finder, LtrDetector) falls back to the bundled binary when it is not on $PATH.
    # Anything still unresolved is left as the bare name and fails loudly at exec time.
    BUNDLED = {
        "ltr_finder": "tools/LTR_FINDER_parallel/bin/LTR_FINDER.x86_64-1.0.7/ltr_finder",
        "LtrDetector": "tools/LtrDetector/bin/LtrDetector",
    }
    for k, v in list(cfg.items()):
        if not v or os.path.isabs(v):
            continue
        if os.sep in v:
            cand = os.path.join(HERE, v)
            if os.path.exists(cand):
                cfg[k] = os.path.abspath(cand)
        elif k in BUNDLED and shutil.which(v) is None:
            cand = os.path.join(HERE, BUNDLED[k])
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
        "--stride", str(a.stride), "--win", str(a.win), "--pad", str(a.pad)], cwd=LIB)
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
              "max_ltr": a.max_ltr, "identity": a.identity,
              # tunables from config.yaml (see load_config); cast here from their strings
              "finder_p": int(cfg.get("ltr_finder_p", 20)),
              "vic": int(cfg.get("ltrharvest_vic", 10)),
              "seed": int(cfg.get("ltrharvest_seed", 20)),
              "trna": (cfg.get("ltr_finder_trna") or None)}

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

    # 4) merge the per-tool .scn files: collapse reciprocally-overlapping calls of the
    #    same locus to a single representative, keeping the highest-priority tool's call.
    #    a.tools order IS the priority (LTR_Finder > LTRharvest > LtrDetector). A call with
    #    no overlap partner forms its own singleton and is always kept, so a locus that only
    #    one tool found (e.g. an LtrDetector-unique element) is never dropped. This both
    #    serves the candidate set (fewer near-duplicates, LTR_Finder-preferred boundaries)
    #    and shrinks the input LTR_retriever adjudicates in the optional intact stage.
    import merge_scn
    merged = os.path.join(d, "candidates.scn")
    present = [s for s in scns if os.path.exists(s)]
    mstats = merge_scn.write_merged(present, merged, reciprocal=a.reciprocal)
    log(f"candidates -> {merged}  ({mstats['out']} candidate LTR-RTs from "
        f"{len(a.tools)} tools; collapsed {mstats['collapsed']} overlap-duplicates "
        f"of {mstats['in']} raw calls, reciprocal>={a.reciprocal})")

    # 4b) candidate deliverables that let ltr_checker replace running the tools separately:
    #     a GFF3 of the merged elements (always) and their sequences as FASTA (unless
    #     --no-candidate-seq). Best-effort: a failure here never discards candidates.scn.
    try:
        import scn_to_fasta_gff as sfg
        gff = os.path.join(d, "candidates.gff3")
        fasta = os.path.join(d, "candidates.fasta")
        want_fa = (not a.no_candidate_seq) and os.path.exists(a.genome + ".fai")
        if (not a.no_candidate_seq) and not want_fa:
            log(f"[candidates] FASTA skipped: no {a.genome}.fai "
                f"(run `samtools faidx {a.genome}`); GFF3 still written")
        nemit = sfg.convert(merged, a.genome, order, fasta, gff, write_fasta=want_fa)
        log(f"candidates -> {gff}" + (f", {fasta}" if want_fa else "")
            + f"  ({nemit} elements)")
    except Exception as ex:   # noqa: BLE001 - deliverables are best-effort
        log(f"WARNING: candidate FASTA/GFF3 emit failed ({ex}); candidates.scn is unaffected")
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
    # lexists (not exists) so a leftover BROKEN symlink from a previous run is detected and
    # replaced; os.symlink would otherwise raise FileExistsError on it.
    if os.path.lexists(g):
        os.remove(g)
    os.symlink(os.path.abspath(a.genome), g)
    # -noanno: skip LTR_retriever's own whole-genome RepeatMasker + LAI pass. It is pure
    # overhead here (~10-20% of the stage) because the separate `annotate` stage runs
    # RepeatMasker with the library anyway; the intact list + LTRlib.fa are still produced.
    # Wrap the call in GNU `/usr/bin/time -v` (when present) so the QF stage records its own
    # wall time and peak RSS -> intact/qf_time.txt, parsed into run_summary.tsv. (On SLURM,
    # sacct MaxRSS captures the whole process tree; time -v captures the LTR_retriever parent.)
    tlog = os.path.join(d, "qf_time.txt")
    # shlex.quote every interpolated path so a space or shell metacharacter in --output
    # (which flows into d/scn/tlog) cannot break or inject into this shell command. The
    # pipe to `tail` is why this one stays shell=True (argv cannot express the pipe).
    tprefix = (f'/usr/bin/time -v -o {shlex.quote(tlog)} '
               if os.path.exists("/usr/bin/time") else "")
    sh(f'cd {shlex.quote(d)} && {tprefix}{shlex.quote(cfg["LTR_retriever"])} '
       f'-genome genome.fa -inharvest {shlex.quote(scn)} '
       f'-noanno -threads {a.threads} 2>&1 | tail -60')
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
    rm = shlex.quote(cfg["RepeatMasker"])
    div = int(cfg.get("rm_div", 40))         # tunable via config.yaml (see load_config)
    cutoff = int(cfg.get("rm_cutoff", 225))
    if a.rm_preset == "edta":
        sh(f'{rm} -e ncbi -pa {pa} -q -no_is -nolow -div {div} '
           f'-lib {shlex.quote(lib_fa)} {shlex.quote(ref)} 2>&1 | tail -30')
    else:  # hite
        sh(f'{rm} -e ncbi -pa {pa} -q -no_is -nolow -div {div} -cutoff {cutoff} -gff '
           f'-lib {shlex.quote(lib_fa)} {shlex.quote(ref)} 2>&1 | tail -30')
    out = ref + ".out"
    if os.path.exists(out):
        log(f"annotate ({a.rm_preset} preset) -> {out}")
    else:
        log("WARNING: RepeatMasker produced no .out")


def _genome_bp(genome):
    """Total genome bp, preferring a .fai index (cheap) over a streamed pass."""
    fai = genome + ".fai"
    if os.path.exists(fai):
        try:
            return sum(int(l.split("\t")[1]) for l in open(fai) if l.strip())
        except (IndexError, ValueError):
            pass
    bp = 0
    with open(genome) as fh:
        for line in fh:
            if not line.startswith(">"):
                bp += len(line.strip())
    return bp


def _count_scn(path):
    if not path or not os.path.exists(path):
        return None
    return sum(1 for l in open(path) if l.strip() and not l.startswith("#"))


def _parse_gnu_time(path):
    """Parse a GNU `/usr/bin/time -v` log -> (wall_s, maxrss_kb), each None if absent."""
    if not path or not os.path.exists(path):
        return None, None
    wall = rss = None
    for l in open(path):
        l = l.strip()
        if l.startswith("Elapsed (wall clock) time"):
            parts = l.split(":", 1)[1].strip().split(":")  # h:mm:ss or m:ss.ss
            try:
                nums = [float(x) for x in parts]
                wall = nums[-1]
                if len(nums) >= 2:
                    wall += nums[-2] * 60
                if len(nums) >= 3:
                    wall += nums[-3] * 3600
            except ValueError:
                pass
        elif l.startswith("Maximum resident set size"):
            try:
                rss = int(l.split(":")[1].strip())
            except (IndexError, ValueError):
                pass
    return wall, rss


def write_run_summary(a, want, timings, regions, scn):
    """R10: write <output>/run_summary.tsv — per-stage wall time, the CNN-flagged fraction,
    and candidates per tool — so users see where time went and whether the gate helped."""
    path = os.path.join(a.output, "run_summary.tsv")
    cand_dir = os.path.join(a.output, "candidates")
    with open(path, "w") as oh:
        oh.write("# ltr_checker run summary\n")
        oh.write(f"genome\t{a.genome}\n")
        oh.write(f"stages\t{','.join(want)}\n")
        oh.write(f"tools\t{','.join(a.tools)}\n")
        oh.write(f"threads\t{a.threads} (finder={a.finder_threads} "
                 f"harvest={a.harvest_threads} detector={a.detector_threads}, "
                 f"{'parallel' if a.parallel_tools else 'serial'})\n")
        # gate: flagged fraction
        if os.path.exists(regions):
            flagged = 0; nreg = 0
            for l in open(regions):
                c = l.split()
                if len(c) >= 3:
                    try:
                        flagged += int(c[2]) - int(c[1]); nreg += 1
                    except ValueError:
                        pass
            try:
                gbp = _genome_bp(a.genome)
                frac = f"{100.0 * flagged / gbp:.1f}%" if gbp else "NA"
            except OSError:
                gbp, frac = 0, "NA"
            oh.write(f"gate_flagged\t{frac}\t({flagged:,} / {gbp:,} bp in {nreg} regions)\n")
        oh.write("#\nstage\twall_s\n")
        for s in want:
            if s in timings:
                oh.write(f"{s}\t{timings[s]:.1f}\n")
        # per-tool candidate counts
        if "candidates" in want:
            oh.write("#\ntool\tcandidates\n")
            for t in a.tools:
                n = _count_scn(os.path.join(cand_dir, f"{t}.scn"))
                oh.write(f"{t}\t{'NA' if n is None else n}\n")
            oh.write(f"merged\t{_count_scn(scn) if _count_scn(scn) is not None else 'NA'}\n")
        # intact count + QF (LTR_retriever) resource use
        if "intact" in want:
            idir = os.path.join(a.output, "intact")
            pl = _retriever_output(idir, ".pass.list")
            ni = _count_scn(pl)
            oh.write(f"#\nintact_LTR_RTs\t{'NA' if ni is None else ni}\n")
            qf_wall, qf_rss = _parse_gnu_time(os.path.join(idir, "qf_time.txt"))
            if qf_wall is not None:
                oh.write(f"qf_wall_s\t{qf_wall:.1f}\n")
            if qf_rss is not None:
                oh.write(f"qf_peak_rss_gb\t{qf_rss / 1048576:.2f}\t"
                         f"({qf_rss:,} kB, LTR_retriever parent process)\n")
    log(f"run summary -> {path}")
    for l in open(path):
        log("  " + l.rstrip())


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
    core.add_argument("--win", type=int, default=50000, choices=[15000, 20000, 25000, 30000, 50000],
                      help="CNN detection window width (bp). 50000 = the bundled model. "
                           "A non-default window needs a gate model trained at that width, "
                           "supplied with --model (not bundled). Smaller windows prune more "
                           "genome on moderate-TE species but cost recall unless --pad compensates.")
    core.add_argument("--pad", type=int, default=-1,
                      help="symmetric safety margin (bp) added to both sides of each candidate "
                           "region before merging. -1 -> auto = max//2 (half the max LTR "
                           "separation), which fully encloses a midpoint-in-window element. "
                           "Use 0 for the pure window union.")
    # structural parameters shared by all three finders
    core.add_argument("--min", type=int, default=1000, help="min LTR separation (element internal), bp")
    core.add_argument("--max", type=int, default=15000, help="max LTR separation (element internal), bp")
    core.add_argument("--min-ltr", type=int, default=100, help="min LTR length, bp")
    core.add_argument("--max-ltr", type=int, default=7000, help="max LTR length, bp")
    core.add_argument("--identity", type=int, default=85, help="min LTR-pair identity (percent)")
    core.add_argument("--reciprocal", type=float, default=0.80,
                      help="reciprocal-overlap fraction at which two tools' calls are treated "
                           "as the same locus and collapsed to the higher-priority tool's call "
                           "when merging the per-tool candidates")
    core.add_argument("--no-candidate-seq", action="store_true",
                      help="skip writing candidates.fasta (candidates.gff3 is still written)")

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
    _default_model = os.path.abspath(os.path.join(HERE, "model", "model_cnn.pth"))
    a.model = os.path.abspath(a.model)
    # The bundled model/model_cnn.pth is the 50 kb gate. A non-default --win needs a gate
    # trained at that window width, which is NOT shipped with the release, so require
    # --model explicitly rather than reaching for an unshipped checkpoint.
    if a.model == _default_model and a.win != 50000:
        sys.exit(f"--win {a.win} needs a gate model trained at that window width, which is "
                 f"not bundled with the release; pass --model /path/to/your_win{a.win}_gate.pth")
    # auto pad = max//2 when not set
    if a.pad < 0:
        a.pad = a.max // 2
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

    timings = {}
    if "gate" in want:
        ts = time.time(); regions = stage_gate(a, cfg); timings["gate"] = time.time() - ts
    if "candidates" in want:
        if not os.path.exists(regions):
            sys.exit("candidates stage needs the gate BED; run --stage gate first")
        ts = time.time(); scn = stage_candidates(a, cfg, regions); timings["candidates"] = time.time() - ts
    if "intact" in want:
        if not os.path.exists(scn):
            sys.exit("intact stage needs candidates.scn; run --stage candidates first")
        ts = time.time(); lib_fa = stage_intact(a, cfg, scn); timings["intact"] = time.time() - ts
    if "annotate" in want:
        if not os.path.exists(lib_fa):
            sys.exit("annotate stage needs the LTR library; run --stage intact first")
        ts = time.time(); stage_annotate(a, cfg, lib_fa); timings["annotate"] = time.time() - ts

    write_run_summary(a, want, timings, regions, scn)
    log(f"done in {time.time() - t0:.1f}s -> {a.output}")


if __name__ == "__main__":
    main()
