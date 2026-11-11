#!/bin/bash
# End-to-end smoke test for ltr_checker on a small bundled genome.
#
# example.fa is a 300 kb slice of rice (Oryza sativa) chromosome 1 (1.35-1.65 Mb)
# that contains real LTR-RTs, so the whole pipeline has something to find:
#
#   gate       CNN flags the LTR-RT-bearing region       (Python + PyTorch only)
#   candidates LTR_Finder + LTRharvest + LtrDetector      (bundled tools + gt)
#   intact     LTR_retriever  -> intact elements + an RM library
#   annotate   RepeatMasker   -> whole-genome LTR-RT annotation
#
# The candidate set (gate + candidates) is the core and needs NO network access.
# The optional downstream stages (intact + annotate) both rely on a FULLY WORKING
# RepeatMasker: LTR_retriever test-runs `RepeatMasker -e ncbi` at startup and aborts
# if it fails, and RepeatMasker itself needs the RMBlast engine + a Dfam library.
# On a sandbox / offline machine where RMBlast or Dfam cannot be fetched, those two
# stages cannot run -- that is an environment limitation, not a problem with
# ltr_checker. So by default this script always runs the core, then attempts the
# downstream stages only if the tools are present and lets them fail gracefully.
#
# Usage:
#   bash test/run_example.sh              # core, then downstream if available (default)
#   bash test/run_example.sh gate         # CNN gate only (Python + PyTorch)
#   bash test/run_example.sh all          # force the full 4-stage run (fail hard if a tool is missing)
#   bash test/run_example.sh gate,candidates   # any explicit --stage subset
set -e
cd "$(dirname "$0")/.."

run() { python ltr_checker.py --genome test/example.fa --output test/out \
                              --name example --threads 4 --device cpu --stage "$1"; }

if [ -n "$1" ]; then
    # explicit stage(s) requested -> run exactly that
    run "$1"
else
    # default: always run the dependency-light core (works everywhere)
    run gate,candidates
    # then try the optional downstream stages, but only if their tools are on PATH,
    # and do not let a RepeatMasker/Dfam failure abort the whole smoke test
    if command -v LTR_retriever >/dev/null 2>&1 && command -v RepeatMasker >/dev/null 2>&1; then
        echo "[run_example] LTR_retriever + RepeatMasker found -> running downstream stages"
        run intact,annotate || echo "[run_example] downstream stages failed (likely RepeatMasker needs RMBlast + a Dfam library) -- the candidate set above is still valid"
    else
        echo "[run_example] LTR_retriever / RepeatMasker not on PATH -> skipping the optional downstream stages (candidate set only)"
    fi
fi

echo
echo "================= results ================="
echo "CNN-flagged regions:"
cat test/out/gate/example_cnn_regions.bed 2>/dev/null || true
if [ -f test/out/candidates/candidates.scn ]; then
    echo "candidate LTR-RTs: $(grep -vc '^#' test/out/candidates/candidates.scn)"
fi
if ls test/out/intact/genome.fa*.pass.list >/dev/null 2>&1; then
    pl=$(ls test/out/intact/genome.fa.mod.pass.list 2>/dev/null || ls test/out/intact/genome.fa.pass.list)
    echo "intact LTR-RTs:    $(grep -vc '^#' "$pl")   ($pl)"
fi
if [ -f test/out/annotate/genome.fa.out ]; then
    echo "genome annotation: test/out/annotate/genome.fa.out"
fi
echo "==========================================="

# Smoke-test assertion: if the candidate stage ran, it MUST produce at least one candidate.
# Without this, a finder that fails silently (empty output) would leave an empty candidates.scn
# and the script would still "pass". example.fa contains real LTR-RTs, so 0 is a real failure.
STAGE="${1:-gate,candidates}"
case ",$STAGE," in
    *,candidates,*|*,all,*)
        scn=test/out/candidates/candidates.scn
        n=0; [ -f "$scn" ] && n=$(grep -vc '^#' "$scn")
        if [ "$n" -lt 1 ]; then
            echo "[run_example] FAIL: candidate stage produced 0 candidates ($scn)" >&2
            exit 1
        fi
        echo "[run_example] OK: $n candidate LTR-RTs"
        ;;
esac
