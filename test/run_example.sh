#!/bin/bash
# End-to-end smoke test for ltr_checker on a small bundled genome.
#
# example.fa is a 300 kb slice of rice (Oryza sativa) chromosome 1 (1.35-1.65 Mb)
# that contains real LTR-RTs, so the WHOLE pipeline has something to find:
#
#   gate       CNN flags the LTR-RT-bearing region
#   candidates LTR_Finder + LTRharvest + LtrDetector  -> a merged candidate set
#   intact     LTR_retriever                          -> intact elements + an RM library
#   annotate   RepeatMasker                           -> whole-genome LTR-RT annotation
#
# The full run takes only a couple of minutes on 4 CPU cores (no GPU needed) and is
# meant as a "does my install work end to end?" check. It needs the external tools
# from environment.yml (gt, LTR_retriever, RepeatMasker) on PATH plus the bundled
# LTR_Finder / LtrDetector under tools/. For a lighter check that needs only Python +
# PyTorch, run just the gate: set STAGE=gate below (or pass it as the first argument).
set -e
cd "$(dirname "$0")/.."

STAGE="${1:-all}"          # default: run everything; pass "gate" for a Python-only check

python ltr_checker.py \
    --genome test/example.fa \
    --output test/out \
    --name example \
    --stage "$STAGE" \
    --threads 4 \
    --device cpu

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
