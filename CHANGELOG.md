# Changelog

## release-cleanup-2026-11 — script-quality pass for a publishable baseline

Correctness, robustness, portability and maintainability cleanup of the release package.
No change to the scientific defaults of the candidate-merge step (owned by a separate
workstream). The end-to-end smoke test (`bash test/run_example.sh gate,candidates`) passes.

### Correctness & robustness
- **run_finders.py** — `ltr_finder` now runs via an argv list (no shell) and its **return
  code is checked**: a missing/incompatible binary (rc 127) now raises instead of being
  silently treated as "no LTR-RTs found" and yielding an empty candidate set. `LtrDetector`
  likewise moved to an argv list (output captured to its log; `check=True` kept).
- **to_scn.py** — the `finder` row parser is now wrapped in `try/except` like the other
  tool branches, so one malformed line can no longer crash the candidates stage (run with
  `check=True`). Dropped malformed rows and elements with inconsistent LTR coordinate order
  are now **counted and reported on stderr** instead of silently discarded.
- **ltr_checker.py** — the `intact` stage no longer crashes on a leftover **broken
  `genome.fa` symlink** (uses `os.path.lexists` + replace).
- **make_input.py** — guards against `--overlap >= --chunk` (which would cut a sequence
  into one piece per base); fails loudly instead.
- **cnn_gate.py** — fixed the standalone `--model` default (pointed at a non-existent
  `lib/model/`; now the package-root `model/`).

### Robustness of shell-outs
- All interpolated paths in the remaining shell commands (`intact` / `annotate` stages)
  are `shlex.quote`-d, so a space or metacharacter in `--output` cannot break or inject.

### Portability (publishable baseline)
- **environment.yml** — added `ltr_finder` (bioconda v1.07) and `samtools` (for the
  `candidates.fasta` `.fai`); removed the unused `pandas` dependency.
- **config.yaml / ltr_checker.py** — tool discovery is now PATH-first with a bundled
  fallback: a bare tool name is taken from `$PATH`, and for the bundled tools
  (`ltr_finder`, `LtrDetector`) falls back to the bundled binary when not on `$PATH`.
  `ltr_finder` is cross-platform via bioconda; only the bundled `LtrDetector` (x86_64)
  remains platform-specific (rebuild script provided).

### Configurability
- **config.yaml** — the previously inlined scientific parameters are now exposed and
  documented: `ltrharvest_vic`, `ltrharvest_seed`, `ltr_finder_p`, `rm_div`, `rm_cutoff`,
  and an optional `ltr_finder_trna` (ltr_finder `-s` PBS tRNA FASTA; off by default, so
  default behaviour is unchanged).
- **ltr_checker.py** — removed the `--win` auto-resolve that pointed at retrain checkpoints
  not shipped in the release; a non-default `--win` now requires an explicit `--model`
  with a clear error.

### Maintainability
- Removed dead code: `model.py` unused encoders/`predict` (+ unused imports),
  `sequence_utils.split_ltr_range` (it still carried the start>end bug it was removed for)
  and a dead empty-case branch, and the unused `look4ltrs` branch in `to_scn.py`.

### Documentation
- **README.md** — documented the portability model (conda tools vs the x86_64 LtrDetector),
  the new config parameters and tRNA option, the `candidates.gff3` / `candidates.fasta`
  outputs, and a **"Large genomes (memory)"** note: the CNN gate loads the whole genome
  into RAM, so very large assemblies should be run per chromosome. The same limit is noted
  at `sequence_utils.genome_dict`.

### Notes
- **GFF3 strand** in `scn_to_fasta_gff.py` was hardcoded `+`; it is now `.` (unknown).
  The `.scn` (LTRharvest `-tabout`) carries no strand column, so `+` was wrong for
  ~half the elements; `.` is the correct value for unknown strand. Full strand propagation
  would require extending the `.scn` contract (deferred; coordinate with the merge redesign).
- The `ltr_finder` tRNA database could not be passed as a single forced default: the bundled
  `tRNAdb/` is a per-species directory and `ltr_finder -s` takes one FASTA, and ltr_finder
  runs cleanly without it. It is therefore exposed as the optional `ltr_finder_trna` knob
  rather than forced on (which would change default candidates and hurt portability).
