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

That is all you need. **`LTR_Finder` and `LtrDetector` ship bundled**
under `tools/` and are used by default — no separate download. `gt`,
`LTR_retriever` and `RepeatMasker` come from the conda environment above.

The bundled LtrDetector is a prebuilt **Linux x86_64** binary. If it will not run on
your platform, rebuild it from the included source:

```bash
bash tools/build_ltrdetector.sh          # needs g++ >= 8.2 and make
```

You normally do not need to touch **`config.yaml`** — its defaults already point at
the bundled tools (relative paths are resolved against the package root, so they
work from any directory). Override a path there (or with an environment variable of
the same name) only to use your own install:

```yaml
ltr_finder:    tools/LTR_FINDER_parallel/bin/LTR_FINDER.x86_64-1.0.7/ltr_finder   # bundled
gt:            gt                                                                 # from conda / $PATH
LtrDetector:   tools/LtrDetector/bin/LtrDetector                                  # bundled
LTR_retriever: LTR_retriever                                                      # from conda / $PATH
RepeatMasker:  RepeatMasker                                                       # from conda / $PATH
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
bash test/run_example.sh         # full pipeline on a small example genome (~2 min, CPU)
bash test/run_example.sh gate    # lighter check: CNN gate only (needs just Python + PyTorch)
```

`test/example.fa` is a 300 kb slice of rice (*Oryza sativa*) chromosome 1 that contains
real LTR-RTs, so the whole chain has something to find. The default run executes **all
four stages** — `gate` → `candidates` → `intact` (LTR_retriever) → `annotate`
(RepeatMasker) — in about two minutes on 4 CPU cores (no GPU needed) and prints a short
summary; on the bundled example it yields 13 candidates, 3 intact LTR-RTs, and a
whole-genome annotation under `test/out/`. It is meant as a "does my install work end to
end?" check, and it needs the external tools from `environment.yml` (`gt`,
`LTR_retriever`, `RepeatMasker`) plus the bundled `LTR_Finder` / `LtrDetector`. The
`gate` argument runs only the CNN gate, which needs nothing but Python + PyTorch.

## Output layout

```
<output>/
├── gate/<name>_cnn_regions.bed        CNN-flagged regions
├── candidates/
│   ├── input_dir/input.fa             region sub-sequences fed to the tools
│   ├── {ltr_finder,ltr_harvest,ltrdetector}.scn
│   └── candidates.scn                 merged candidates (LTR_retriever input)
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

## Requirements

- Python ≥ 3.9, PyTorch, NumPy, Biopython, pandas, tqdm
- From conda (`environment.yml`): genometools (`gt`), LTR_retriever, RepeatMasker, perl
- Bundled under `tools/` (no download): LTR_Finder, LtrDetector

## Citation

If you use ltr_checker, please cite it (see `CITATION.cff`).

## License

MIT — see `LICENSE`.
