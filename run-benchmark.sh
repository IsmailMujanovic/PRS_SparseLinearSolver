#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

if ! command -v make >/dev/null 2>&1; then
    echo "error: make is required" >&2
    exit 1
fi

if ! command -v c++ >/dev/null 2>&1; then
    echo "error: a C++ compiler is required" >&2
    exit 1
fi

SUPERLU_AVAILABLE=0

if command -v brew >/dev/null 2>&1; then
    SUPERLU_PREFIX="$(brew --prefix superlu 2>/dev/null || true)"
    OPENBLAS_PREFIX="$(brew --prefix openblas 2>/dev/null || true)"
    if [ -n "$SUPERLU_PREFIX" ] && [ -f "$SUPERLU_PREFIX/include/slu_zdefs.h" ] && [ -n "$OPENBLAS_PREFIX" ]; then
        SUPERLU_AVAILABLE=1
    fi
fi

if [ "$SUPERLU_AVAILABLE" -eq 0 ]; then
    if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists superlu; then
        SUPERLU_AVAILABLE=1
    fi
fi

if [ "$SUPERLU_AVAILABLE" -eq 0 ]; then
    if [ -f /usr/include/slu_zdefs.h ] || [ -f /usr/local/include/slu_zdefs.h ]; then
        SUPERLU_AVAILABLE=1
    fi
fi

if [ "$SUPERLU_AVAILABLE" -eq 0 ]; then
    if [ -f /usr/include/superlu/slu_zdefs.h ] || [ -f /usr/local/include/superlu/slu_zdefs.h ]; then
        SUPERLU_AVAILABLE=1
    fi
fi

if [ "$SUPERLU_AVAILABLE" -eq 0 ]; then
    echo "error: SuperLU is not available on this machine." >&2
    echo "This script does not download dependencies. Install SuperLU once, then rerun." >&2
    exit 1
fi

make -s bench

if [ "$#" -eq 0 ]; then
    set -- --repeat 3 --large --grid 50 --grid 60 --grid 80 --grid 100 --grid 120
fi

exec ./benchmark/run.sh "$@"
