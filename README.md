# ltr_checker

**ltr_checker builds an LTR-retrotransposon (LTR-RT) candidate set for a plant genome
with a CNN gate + three structural tools.** A convolutional neural network scans the
genome in 50 kb windows and flags the regions likely to contain LTR-RTs; three structural
tools (**LTR_Finder · LTRharvest · LtrDetector**) then run **only inside those
regions**, and their calls are remapped to genome coordinates and merged into one
candidate set. That candidate set is the tool's primary output and what it produces by
default. (LTR_Finder is driven directly, with its own CLI parameters; the piece-splitting
parallelization follows the LTR_FINDER_parallel strategy.)

Two **optional downstream stages** are kept for convenience and enabled with `--stage`:
`LTR_retriever` turns the candidate set into intact elements + a RepeatMasker library, and
`RepeatMasker` annotates the whole genome with that library.

```
              ┌─────────────┐   50 kb sliding windows
   genome ───►│  CNN gate   │──────────────────────────►  flagged regions (BED)
              └─────────────┘
                     │  run only inside flagged regions
                     ▼
   ┌────────────────────────────────────────────────┐
   │ LTR_Finder · LTRharvest · LtrDetector           │ ──► candidate LTR-RTs (candidates.scn)
   └────────────────────────────────────────────────┘        ◄── default output

   ─ ─ ─ optional downstream stages (--stage intact,annotate | all) ─ ─ ─
                     │
                     ▼
              ┌───────────────┐
              │ LTR_retriever │ ──►  intact LTR-RTs (pass.list)  +  RM library (LTRlib.fa)
              └───────────────┘
                     │
                     ▼
              ┌───────────────┐
              │ RepeatMasker  │ ──►  whole-genome LTR-RT annotation (.out / .gff)
              │ EDTA | HiTE   │
              └───────────────┘
```

## Stages

| stage        | what it does                                                        | key output | default |
|--------------|---------------------------------------------------------------------|------------|---------|
| `gate`       | CNN flags LTR-RT-bearing regions (50 kb windows, 10 kb stride)      | `gate/<name>_cnn_regions.bed` | ✓ |
| `candidates` | LTR_Finder + LTRharvest + LtrDetector inside those regions          | `candidates/candidates.scn` | ✓ |
| `intact`     | *(optional)* LTR_retriever → intact elements + RepeatMasker library | `intact/genome.fa.mod.pass.list`, `intact/genome.fa.mod.LTRlib.fa` | |
| `annotate`   | *(optional)* RepeatMasker over the genome with that library         | `annotate/genome.fa.out` | |

**The default is `--stage gate,candidates`** (the candidate set). Add the downstream
stages with `--stage intact,annotate`, or run everything with `--stage all`.

## Installation

```bash
git clone https://github.com/morningsun77/ltr_checker.git
cd ltr_checker
conda env create -f environment.yml
conda activate ltr_checker
```

That is all you need. The conda environment provides `ltr_finder` (bioconda), `gt`,
`samtools`, and the optional `LTR_retriever` / `RepeatMasker`. **`LtrDetector` ships
bundled** under `tools/` (it is not on bioconda) and is used automatically. A prebuilt
`ltr_finder` is also bundled and used as a fallback if it is not on `PATH`.

The bundled LtrDetector is a prebuilt **Linux x86_64** binary. If it will not run on
your platform, rebuild it from the included source:

```bash
bash tools/build_ltrdetector.sh          # needs g++ >= 8.2 and make
```

**Platform note.** The conda tools (`ltr_finder`, `gt`, `samtools`, …) install on any
platform bioconda supports. The only non-portable piece is the bundled LtrDetector binary
(Linux x86_64), which has the rebuild script above for other platforms.

You normally do not need to touch **`config.yaml`**. Each tool path is resolved as: an
absolute path as-is; a relative path (with a `/`) against the package root; a bare binary
name from `PATH`, and for the bundled tools (`ltr_finder`, `LtrDetector`) falling back to
the bundled binary when the name is not on `PATH`. An environment variable of the same name
wins over the file. `config.yaml` also exposes the per-tool scientific parameters:

```yaml
# tool paths (bare name = $PATH, with bundled fallback for ltr_finder / LtrDetector)
ltr_finder:    ltr_finder        # bioconda, or bundled x86_64 fallback
gt:            gt                 # from conda / $PATH
LtrDetector:   LtrDetector        # bundled under tools/ (not on bioconda)
LTR_retriever: LTR_retriever      # from conda / $PATH
RepeatMasker:  RepeatMasker       # from conda / $PATH
# tunable parameters (otherwise hidden in the tool invocations)
ltrharvest_vic:    10             # LTRharvest -vic / -seed
ltrharvest_seed:   20
ltr_finder_p:      20             # ltr_finder -p
ltr_finder_trna:                  # optional ltr_finder -s tRNA FASTA (empty = off); e.g.
                                  #   tools/LTR_FINDER_parallel/bin/.../tRNAdb/Os-tRNAs.fa
rm_div:            40             # RepeatMasker -div / -cutoff (annotate stage)
rm_cutoff:        225
```

## Usage

```bash
# default: CNN + three tools -> candidate set (candidates.scn)
python ltr_checker.py --genome genome.fa --output result --threads 16

# use the GPU for the CNN gate
python ltr_checker.py --genome genome.fa --output result --device cuda

# run the three tools concurrently instead of one after another
python ltr_checker.py --genome genome.fa --output result --parallel-tools

# also run the optional downstream stages (LTR_retriever + RepeatMasker)
python ltr_checker.py --genome genome.fa --output result --stage all
python ltr_checker.py --genome genome.fa --output result --stage all --rm-preset hite
```

### Candidate-set options (the core: CNN gate + three tools)

| option        | default                              | meaning |
|---------------|--------------------------------------|---------|
| `--genome`    | *(required)*                         | genome FASTA |
| `--output`    | *(required)*                         | output directory |
| `--stage`     | `gate,candidates`                    | stages to run: comma subset of `gate,candidates,intact,annotate`, or `all` |
| `--tools`     | `ltr_finder,ltr_harvest,ltrdetector` | candidate tools to run |
| `--threads`   | `10`                                 | unified thread count for every tool |
| `--parallel-tools` | *(off)*                         | run the three tools concurrently (peak threads = sum) instead of one after another |
| `--finder-threads` | `--threads`                     | parallel LTR_Finder workers |
| `--harvest-threads` | `--threads`                    | parallel LTRharvest workers |
| `--detector-threads` | `--threads`                   | threads for LtrDetector |
| `--device`    | `cpu`                                | `cpu` or `cuda` (CNN gate) |
| `--model`     | `model/model_cnn.pth`                | CNN weights |
| `--stride`    | `10000`                              | CNN window stride (bp) |

Structural parameters (**honored by all three tools** — LTR_Finder too):

| option        | default   | meaning |
|---------------|-----------|---------|
| `--min` / `--max`       | `1000` / `15000` | LTR separation (element internal span), bp |
| `--min-ltr` / `--max-ltr` | `100` / `7000` | LTR length, bp |
| `--identity`  | `85`      | min LTR-pair identity (percent) |

### Optional downstream stages (LTR_retriever / RepeatMasker)

Enabled only when `intact` / `annotate` are requested via `--stage`.

| option        | default   | meaning |
|---------------|-----------|---------|
| `--rm-preset` | `edta`    | RepeatMasker preset for the `annotate` stage: `edta` or `hite` |
| `--config`    | `config.yaml` | external-tool paths (also used by the tools) |

### Threading

Every tool is driven by `--threads` (default 10). The `candidates` stage runs the three
tools **one after another by default**, each using all `--threads`; pass
`--parallel-tools` to run them **concurrently** (then peak concurrency is the sum of the
three per-tool counts, which default to `--threads` each — set `--finder-threads` /
`--harvest-threads` / `--detector-threads` to tune). Serial keeps peak memory to one tool
at a time (recommended for large genomes); parallel trades memory for wall-clock. Under
the hood the flagged regions are cut into overlapping 5 Mb pieces and the pieces are
balanced into per-tool worker bundles: LTR_Finder and LTRharvest each run one process per
bundle (one `gt` index per bundle for LTRharvest), and LtrDetector uses its own
`-nThreads`. The piece-splitting mirrors the LTR_FINDER_parallel strategy. Downstream,
`LTR_retriever` and `RepeatMasker` also use `--threads`.

Stages are resumable: each stage reads the previous stage's output from `--output`,
so you can rerun `--stage annotate` alone after `intact` has finished.

### RepeatMasker presets

| preset | command |
|--------|---------|
| `edta` | `RepeatMasker -e ncbi -pa N -q -no_is -nolow -div 40 -lib LTRlib.fa genome.fa` |
| `hite` | `RepeatMasker -e ncbi -pa N -q -no_is -nolow -div 40 -cutoff 225 -gff -lib LTRlib.fa genome.fa` |

## Quick test

```bash
bash test/run_example.sh         # core, then downstream if available (~2 min, CPU)
bash test/run_example.sh gate    # lighter check: CNN gate only (needs just Python + PyTorch)
bash test/run_example.sh all     # force the full 4-stage run (fails hard if a tool is missing)
```

`test/example.fa` is a 300 kb slice of rice (*Oryza sativa*) chromosome 1 that contains
real LTR-RTs, so the whole chain has something to find. By default the script first runs
the **candidate set** (`gate` → `candidates`) — which needs no network access — and then
attempts the optional downstream stages (`intact` → `annotate`) only if `LTR_retriever`
and `RepeatMasker` are on `PATH`, letting them fail gracefully. On a complete install it
finishes in about two minutes on 4 CPU cores (no GPU needed) and yields a handful of
candidate LTR-RTs (the script asserts at least one), a few intact LTR-RTs, and a
whole-genome annotation under `test/out/`. The `gate` argument runs only the CNN gate
(Python + PyTorch, nothing else); `all` forces every stage.

> **About the downstream stages and RepeatMasker.** `LTR_retriever` test-runs
> `RepeatMasker -e ncbi` at startup and aborts if it fails, so `intact`/`annotate` need a
> working RepeatMasker with the **RMBlast** engine. You do **not** download these by hand:
> the `conda env create -f environment.yml` above installs `repeatmasker`, which pulls in
> **RMBlast automatically** (it is a package dependency) and ships with a **minimal Dfam
> library** that is enough for RepeatMasker to run. ltr_checker's `annotate` stage masks
> the genome with the **custom `-lib` built from your own candidates** (`LTRlib.fa`), so it
> does **not** require the full Dfam download — a complete Dfam is optional and only
> matters if you separately want genome-wide homology masking against all known TE
> families. The stages only truly fail on an **offline/sandboxed** machine where conda
> cannot install RMBlast at all (e.g. `repeatmasker.org`, NCBI and `dfam.org` blocked);
> that is an environment limitation, not a problem with ltr_checker, and the candidate set
> still runs fine. Do **not** substitute the experimental Rust RMBlast port from GitHub;
> its authors warn against pairing it with current RepeatMasker.

## Output layout

```
<output>/
├── run_summary.tsv                    per-stage wall time, CNN-flagged fraction, candidates per tool
├── gate/<name>_cnn_regions.bed        CNN-flagged regions
├── candidates/
│   ├── input_dir/input.fa             region sub-sequences fed to the tools
│   ├── {ltr_finder,ltr_harvest,ltrdetector}.scn
│   ├── candidates.scn                 merged candidates (LTR_retriever input)
│   ├── candidates.gff3                merged candidates as annotations (always written)
│   └── candidates.fasta               merged candidate sequences (only if <genome>.fai exists;
│                                      run `samtools faidx <genome>` first, or --no-candidate-seq)
│
│   ── the two directories below are written only when the optional stages run ──
├── intact/
│   ├── genome.fa.mod.pass.list        intact LTR-RTs (LTR_retriever adds the `.mod` infix)
│   └── genome.fa.mod.LTRlib.fa        RepeatMasker library
└── annotate/genome.fa.out             genome-wide LTR-RT annotation
```

## Model & training data

The bundled CNN weights are `model/model_cnn.pth`. The gate is a convolutional network
trained on 50 kb windows (a positive window carries one LTR-RT padded with mRNA; a
negative window is mRNA only), over a balanced set of 241,171 positive and 241,171
negative instances drawn from full-length LTR-RTs across broad plant phylogenetic
coverage (angiosperms, gymnosperms, and others).

The **LTR-RTs (and mRNAs) used to train and test the model are publicly available
on Zenodo**:

> https://zenodo.org/records/10893270

## Files

```
ltr_checker.py       pipeline driver (this is the entry point)
config.yaml          external-tool paths (defaults point at tools/)
model/model_cnn.pth  trained CNN weights
lib/                 supporting modules (CNN gate, coordinate remapping, I/O)
tools/               bundled LTR_Finder (LTR_FINDER_parallel pkg) + LtrDetector (+ build_ltrdetector.sh)
environment.yml      conda environment
test/                small example genome + run_example.sh
```

## Large genomes (memory)

The CNN gate loads the **whole genome into memory** (as Python strings), so peak RAM scales
with genome size (roughly 1 byte/bp plus overhead). A multi-Gb plant genome can exhaust
memory. For such assemblies, run ltr_checker **per chromosome** (split the FASTA and run
each sequence separately), then concatenate the resulting `candidates.scn`. The optional
`intact`/`annotate` stages (LTR_retriever / RepeatMasker) are not affected by this.

## Requirements

- Python ≥ 3.9, PyTorch, NumPy, Biopython, tqdm
- From conda (`environment.yml`): `ltr_finder`, genometools (`gt`), `samtools`, and the
  optional-stage tools LTR_retriever, RepeatMasker, perl
- Bundled under `tools/` (no download): LtrDetector (Linux x86_64; rebuild script included).
  A prebuilt `ltr_finder` is also bundled as a fallback when it is not on `PATH`.
- Only for the optional `intact`/`annotate` stages: a working RepeatMasker with the
  **RMBlast** engine. Both come from the conda install — `repeatmasker` pulls in RMBlast
  automatically and bundles a minimal Dfam, and the `annotate` stage uses your own
  `LTRlib.fa` as the library, so no separate RMBlast/Dfam download is needed (only a
  network connection for conda itself). The candidate set (`gate`/`candidates`) needs none
  of this.

## Citation

If you use ltr_checker, please cite it (see `CITATION.cff`).

## License

MIT — see `LICENSE`.
