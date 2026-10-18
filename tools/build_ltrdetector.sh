#!/bin/bash
# Rebuild the bundled LtrDetector binary from source.
# Only needed if the prebuilt tools/LtrDetector/bin/LtrDetector does not run on
# your platform (it is a Linux x86_64 ELF). Requires GNU g++ >= 8.2 and make.
#   conda install -c conda-forge make cxx-compiler    # if you lack a compiler
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$HERE/LtrDetector/src"
CXX="${CXX:-g++}"

echo "[build_ltrdetector] compiler: $CXX ($($CXX --version 2>/dev/null | head -1))"
cd "$SRC"
make bin
make tr -j CXX="$CXX"
echo "[build_ltrdetector] done -> $HERE/LtrDetector/bin/LtrDetector"
"$HERE/LtrDetector/bin/LtrDetector" 2>&1 | head -1 || true
