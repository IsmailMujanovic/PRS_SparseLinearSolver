# Benchmark: serial complex sparse Cholesky solver vs SuperLU 7.0.1

This document records the complete benchmark of the project's four-phase serial complex sparse Cholesky solver against SuperLU. It contains the methodology, the exact phase mapping, all measured results, the bottleneck analysis, and the reproduction instructions.

The benchmark was run on **2026-11-09**.

## 1. Goal

The solver in this repository follows the four-step direct solution described in the course material:

1. **Ordering** — find a permutation `P` that reduces fill in the Cholesky factor.
2. **Symbolic factorization** — determine the structure of `L` and allocate storage.
3. **Numeric factorization** — compute `PAP^H = L L^H`.
4. **Triangular solution** — solve `L y = P b`, `L^H z = y`, and `x = P^H z`.

The benchmark measures each phase separately and compares it with the closest corresponding phase in SuperLU. The goal is not only to report total speed, but to identify which phase is responsible for the gap so that optimization effort can be directed correctly.

SuperLU is used as the reference because it is a mature, highly optimized sparse direct solver. It is not a mathematically optimal solver: it computes an LU factorization of a general unsymmetric matrix, while this project computes a Cholesky factorization of a Hermitian positive-definite matrix. The comparison therefore accounts for the algorithmic difference between LU and Cholesky.

### Key finding

The Minimum Degree ordering implementation is the bottleneck. At the largest measured problem (`n = 14400`), ordering consumes about 94% of the project's runtime and essentially all of the gap to SuperLU. Symbolic factorization and triangular solution together are below 1.5% of the runtime. Numeric Cholesky factorization is competitive: with the same permutation it is faster in wall-clock time than SuperLU's LU factorization, despite SuperLU having a somewhat higher flop rate. The detailed evidence is in section 8.

## 2. Environment

| Item | Value |
|---|---|
| Machine | Apple M1, 8 GiB RAM |
| Operating system | macOS 27.2, build 26B123, arm64 |
| Compiler | Apple clang 21.0.0 (clang-2100.4.1.1), target arm64-apple-darwin27.2.0 |
| Compiler flags | `-std=c++17 -O3 -Wall -Wextra -Iinclude` |
| SuperLU | 7.0.1, installed with Homebrew |
| BLAS | OpenBLAS 0.3.34, installed as a SuperLU dependency |
| Threading | `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1` |
| Timing | wall-clock, median of 3 measured repetitions after one SuperLU warm-up |
| Build type | optimized, no debug assertions |

SuperLU was installed with:

```bash
brew install superlu
```

The benchmark executable and the project solver use the same compiler, optimization level, and complex type (`std::complex<double>` in the project, `doublecomplex` in SuperLU).

All timings are machine-specific. Absolute times will differ on other CPUs, compilers, BLAS builds, and thread settings.

## 3. Reproduction

The default benchmark output is intentionally small: one table per matrix with only the four phase times and the total, comparing the project solver with SuperLU.

```bash
./run-benchmark.sh
```

The script checks that SuperLU is available, builds the benchmark if necessary, and runs the full suite with `--repeat 3`. It does not download dependencies; SuperLU and OpenBLAS must already be installed.

The detailed tables in this document were produced with the diagnostic mode enabled on the full suite:

```bash
./run-benchmark.sh --repeat 3 --large \
  --grid 50 --grid 60 --grid 80 --grid 100 --grid 120 \
  --detailed --csv build/benchmark_results.csv
```

On Windows, use `run-benchmark.bat` or `run-benchmark.ps1` from PowerShell. The complete per-platform instructions are in [running.md](running.md).

Manual build and run:

```bash
make program
make bench
```

Verify that the project solver still passes its own test suite:

```bash
./build/program
```

The last line must be `ALL OK`.

Run the default benchmark suite:

```bash
./benchmark/run.sh --repeat 3
```

Run the complete suite used for this document, including the larger 2D problems and the diagnostic tables:

```bash
./benchmark/run.sh \
  --repeat 3 \
  --large \
  --grid 50 \
  --grid 60 \
  --grid 80 \
  --grid 100 \
  --grid 120 \
  --detailed \
  --csv build/benchmark_results.csv
```

Available command-line options:

| Option | Meaning |
|---|---|
| `--repeat N` | Number of measured repetitions; the median is reported. Default: 3. |
| `--large` | Add `random-900` and `laplacian2d-40` to the default suite. |
| `--grid K` | Add a complex 2D Laplacian on a `K x K` grid. |
| `--random N:DENSITY:SEED` | Add a random HPD matrix. |
| `--detailed` | Print the extra diagnostic tables and accuracy information. |
| `--csv FILE` | Write one row per matrix and phase to a CSV file. |
| `--help` | Print the usage message. |

The CSV format is:

```text
matrix,n,phase,ours_ms,superlu_ms,ours_minus_superlu_ms,ours_over_superlu_x
```

The difference is in milliseconds: positive means the project phase is slower, negative means it is faster. The ratio column is above 1 when SuperLU is faster and below 1 when the project phase is faster.

The benchmark script forces single-threaded BLAS so that the comparison is not affected by OpenBLAS thread scheduling.

On Linux, if SuperLU and OpenBLAS are installed system-wide, override the detected prefixes:

```bash
make bench SUPERLU_PREFIX=/usr OPENBLAS_PREFIX=/usr
```

## 4. Phase definitions and mapping to SuperLU

The project phases are timed directly in `benchmark/bench_superlu.cpp`. SuperLU does not expose all four phases through a single public timer, so the mapping is:

| Phase | Project implementation | SuperLU reference | Measurement |
|---|---|---|---|
| 1. Ordering | `MD(A)` | `get_perm_c(MMD_AT_PLUS_A)` and `get_perm_c(COLAMD)` | Direct calls, timed separately |
| 2. Symbolic | `permute(A, P)` and `symbolicFactorization(PAPt)` | `ETREE` plus estimated symbolic work | `ETREE + (FACT_first - FACT_reuse)` |
| 3. Numeric | `numericFactorization(PAPt, symbolic)` | numeric-only `FACT` | Second `zgssvx` call with `Fact = SamePattern_SameRowPerm` |
| 4. Triangular solution | `solve(L, P, b)` | `SOLVE` | `zgstrs` timing from the first `zgssvx` call |

### 4.1 Why phase 2 is estimated

SuperLU performs the column symbolic factorization inside `zgstrf`; it does not expose a separate public symbolic-factorization timer. The only public analysis timer is `ETREE`, which measures `sp_preorder` (forming the permuted matrix and the elimination tree) but not the subsequent column DFS and storage setup.

To obtain a phase-aligned estimate:

1. The first `zgssvx` call uses `Fact = DOFACT`; its `FACT` time contains SuperLU's symbolic and numeric work.
2. A second `zgssvx` call uses `Fact = SamePattern_SameRowPerm`. It reuses the column permutation, row permutation, elimination tree, and the allocated `L`/`U` structures, so it performs essentially numeric work only.
3. `FACT_reuse` is used as the numeric factor time.
4. The symbolic factor time is estimated as:

```text
SuperLU symbolic = ETREE + max(0, FACT_first - FACT_reuse)
SuperLU numeric  = min(FACT_reuse, FACT_first)
SuperLU total    = COLPERM + SuperLU symbolic + SuperLU numeric + SOLVE
```

This split is valid only when row pivoting is disabled, because `SamePattern_SameRowPerm` also reuses the row permutation. The phase-aligned SuperLU reference therefore uses:

| Option | Value |
|---|---|
| `Equil` | `NO` |
| `IterRefine` | `NOREFINE` |
| `DiagPivotThresh` | `0.0` |
| `SymmetricMode` | `YES` |
| `ColPerm` | `MMD_AT_PLUS_A` |

With no row pivoting, `L` and `U` are structurally symmetric for a Hermitian matrix, and the second factorization is a clean numeric-only measurement. The estimate still includes warm-cache and preallocated-memory effects, so it is an approximation, not an independent symbolic timer.

SuperLU's stock configuration is also measured for context:

| Option | Value |
|---|---|
| `Equil` | `NO` |
| `IterRefine` | `NOREFINE` |
| `DiagPivotThresh` | `1.0` (partial pivoting) |
| `ColPerm` | `COLAMD` |

The stock configuration is not used for the phase ratios because partial pivoting and COLAMD produce substantially more fill on these HPD matrices. It is reported as the out-of-the-box comparison.

A third SuperLU mode, **same-order**, feeds the project's MD permutation to SuperLU through `ColPerm = MY_PERMC`, `DiagPivotThresh = 0.0`, and `SymmetricMode = YES`. This mode removes ordering differences and isolates the numeric kernel comparison. Its `COLPERM` time is zero by construction because the ordering is supplied by the caller.

### 4.2 Treatment of `permute`

`permute(A, P)` forms the lower triangle of `PAP^H`. It is counted as part of phase 2 because it prepares the matrix structure used by `symbolicFactorization`. Its time is also reported separately as `permute=` in the benchmark output. At the largest case (`n = 14400`), `permute` takes 1.426 ms out of 2.954 ms of phase 2, and is therefore a meaningful part of the symbolic phase.

## 5. Test matrices

All benchmark matrices are complex Hermitian positive-definite by construction. The benchmark does not call `validateHPD` because that validation performs a dense `O(n^3)` Cholesky check, which would dominate the measurement for large `n`; the regular test program `./build/program` does use it.

| Name | `n` | `nnz(A)` | Construction |
|---|---:|---:|---|
| `tridiagonal-300` | 300 | 898 | complex tridiagonal, diagonal 4, subdiagonal `1 - i`, superdiagonal `1 + i` |
| `random-60` | 60 | 350 | random HPD, density 0.08, seed 1 |
| `random-300` | 300 | 2198 | random HPD, density 0.02, seed 2 |
| `laplacian2d-20` | 400 | 1920 | complex 2D Laplacian, 20 x 20 grid |
| `laplacian2d-30` | 900 | 4380 | complex 2D Laplacian, 30 x 30 grid |
| `random-900` | 900 | 2628 | random HPD, density 2/900, seed 3 |
| `laplacian2d-40` | 1600 | 7840 | complex 2D Laplacian, 40 x 40 grid |
| `custom-grid-50` | 2500 | 12300 | complex 2D Laplacian, 50 x 50 grid |
| `custom-grid-60` | 3600 | 17760 | complex 2D Laplacian, 60 x 60 grid |
| `custom-grid-80` | 6400 | 31680 | complex 2D Laplacian, 80 x 80 grid |
| `custom-grid-100` | 10000 | 49600 | complex 2D Laplacian, 100 x 100 grid |
| `custom-grid-120` | 14400 | 71520 | complex 2D Laplacian, 120 x 120 grid |

### 5.1 Random HPD generator

For each lower-triangular position `(i, j)` with `j < i`, the generator draws a uniform number in `[0, 1)`. If it is below the requested density, it sets

```text
a[i][j] = x + i y,      x, y ~ Uniform(-1, 1)
a[j][i] = conj(a[i][j])
```

and leaves the entry zero otherwise. The diagonal is set to

```text
a[i][i] = 1 + sum_j |a[i][j]|.
```

The matrix is strictly diagonally dominant with a positive real diagonal, hence Hermitian positive-definite. The full Hermitian matrix is stored in CSC form, with both triangles present. This matches the project's existing `randomHPD` generator.

### 5.2 Complex 2D Laplacian generator

For a `k x k` grid, the matrix has `n = k^2`. Each grid edge carries a complex unit-magnitude coupling. For an undirected edge between `a` and `b`, one orientation stores

```text
weight = -exp(i theta)
```

and the opposite orientation stores its conjugate. The phases are:

```text
horizontal: theta = 0.37 + 0.05 * ((i + j) mod 5)
vertical:   theta = -0.81 + 0.07 * ((i + 2 j) mod 5)
```

The diagonal is `1 + sum |a[i][j]|` over all incident edges. The resulting matrix is Hermitian, strictly diagonally dominant, and genuinely complex. It has the usual five-point Laplacian sparsity pattern with a complex phase on each off-diagonal.

### 5.3 Right-hand side and known solution

For every matrix the benchmark uses the known solution

```text
x_known[i] = (i + 1) + i (0.5 i - 1)
```

and computes `b = A x_known`. The same `A` and `b` are passed to the project solver and to every SuperLU mode.

## 6. Measurement methodology

### 6.1 Repetitions

Each solver configuration is run `repeat = 3` times. The reported phase time is the median of the three runs. SuperLU receives one additional unmeasured warm-up run per configuration to remove first-call library and allocator overhead. The project solver does not receive a separate warm-up because it has no external library initialization; its three runs are measured directly.

The median reduces the influence of scheduler noise, page faults, and allocator behavior. For sub-millisecond cases, timer resolution and process overhead still dominate; the large 2D cases are the meaningful performance measurements.

Absolute timings are inherently noisy. The numbers in this document come from one full benchmark run; re-running the same command on the same machine typically changes absolute times by several percent, especially for sub-millisecond phases and for the SuperLU run that is warmed up only once per configuration. The qualitative conclusions are based on the scaling trend across many sizes, not on a single timing sample.

### 6.2 Accuracy

For a computed solution `x`, the benchmark reports:

```text
relative residual = ||b - A x||_inf / (||A||_inf ||x||_inf + ||b||_inf)
forward error     = ||x - x_known||_inf
```

The relative residual measures backward stability. The forward error also accounts for the conditioning of `A`. Every solver run must produce `SuperLU info = 0` or the warning `n + 1`; the warning is recorded but not treated as a failure. The benchmark itself fails if the relative residual exceeds `1e-8`.

### 6.3 Fill, factor storage, and work

| Quantity | Project solver | SuperLU |
|---|---|---|
| Factor nonzeros | `nnz(L)` including diagonal | `nnz(L) + nnz(U) - n`, counting the diagonal once |
| Cholesky-equivalent nonzeros | `nnz(L)` | `(nnz(L) + nnz(U)) / 2` in no-pivot mode |
| Fill ratio | `nnz(L) / nnz(A)` | `(nnz(L) + nnz(U) - n) / nnz(A)` |
| Factor storage | `nnz(L) * 16 + n * 8` bytes | `mem_usage.for_lu` bytes |
| Numeric work | estimated real flops from the cmod pattern | `stat->ops[FACT]` |

The Cholesky-equivalent conversion is exact for the no-pivot HPD mode when `L` and `U` are structurally symmetric: `nnz(L) + nnz(U) = 2 (nnz(L_chol) - n) + 2 n`, so `(nnz(L) + nnz(U)) / 2 = nnz(L_chol)`. This is the fair way to compare SuperLU's fill with the project's Cholesky fill.

The project's numeric flop estimate is:

```text
flops ≈ 8 * (number of complex multiply-add updates in cmod)
        + 3 * (number of off-diagonal nonzeros of L)
        + 6 * n
```

The factor 8 corresponds to four real multiplications and four real additions per complex multiply-add. The estimate is approximate and excludes some lower-order operations. SuperLU's count is the solver's own reported operation count. GFLOP/s is computed as `flops / (numeric_seconds * 1e9)`.

## 7. Results

### 7.1 Headline case: complex 2D Laplacian, n = 14400

This is the largest measured case. SuperLU uses the phase-aligned HPD mode (MMD ordering, no row pivoting).

| Phase | Project (ms) | SuperLU (ms) | Slowdown | Project share | Gap (ms) | Gap share |
|---|---:|---:|---:|---:|---:|---:|
| 1. Ordering | 401.197 | 1.867 | 214.9x | 94.3% | 399.330 | 100.3% |
| 2. Symbolic (+ permute) | 2.954 | 0.998 | 3.0x | 0.7% | 1.955 | 0.5% |
| 3. Numeric | 19.114 | 23.816 | 0.80x | 4.5% | -4.701 | -1.2% |
| 4. Triangular solve | 1.626 | 0.844 | 1.9x | 0.4% | 0.782 | 0.2% |
| **Total** | **425.479** | **27.416** | **15.5x** | **100.0%** | **398.063** | **100.0%** |

Additional facts for this case:

| Metric | Project | SuperLU HPD |
|---|---:|---:|
| `nnz(A)` | 71520 | 71520 |
| Factor nonzeros | 342036 (Cholesky `L`) | 563958 (`L + U - n`) |
| Cholesky-equivalent nonzeros | 342036 | 289179 |
| Factor storage | 5.109 MiB | 10.151 MiB |
| Numeric flops | 103.9 million | 143.6 million |
| Numeric rate | 5.44 GFLOP/s | 6.03 GFLOP/s |
| Relative residual | 7.77e-16 | 5.08e-16 |
| Forward error | 3.71e-11 | 2.27e-11 |

The ordering phase is 94.3% of the project's runtime and accounts for essentially the entire gap. The numeric phase is actually faster in wall-clock time than SuperLU's numeric phase, despite a lower flop rate, because Cholesky performs less arithmetic than LU.

### 7.2 Second headline case: complex 2D Laplacian, n = 10000

| Phase | Project (ms) | SuperLU (ms) | Slowdown | Project share | Gap (ms) | Gap share |
|---|---:|---:|---:|---:|---:|---:|
| 1. Ordering | 215.973 | 1.285 | 168.1x | 93.7% | 214.689 | 100.0% |
| 2. Symbolic (+ permute) | 1.976 | 0.453 | 4.4x | 0.9% | 1.523 | 0.7% |
| 3. Numeric | 11.573 | 13.077 | 0.89x | 5.0% | -1.504 | -0.7% |
| 4. Triangular solve | 1.085 | 0.494 | 2.2x | 0.5% | 0.591 | 0.3% |
| **Total** | **230.607** | **15.963** | **14.4x** | **100.0%** | **214.644** | **100.0%** |

| Metric | Project | SuperLU HPD |
|---|---:|---:|
| Factor nonzeros | 221032 | 361346 (`L + U - n`) |
| Cholesky-equivalent nonzeros | 221032 | 185673 |
| Factor storage | 3.296 MiB | 6.533 MiB |
| Numeric flops | 59.4 million | 77.6 million |
| Numeric rate | 5.14 GFLOP/s | 5.93 GFLOP/s |
| Relative residual | 8.40e-16 | 4.79e-16 |

The pattern is the same as at `n = 14400`: ordering owns 93.7% of the runtime and 100% of the gap, while numeric factorization is faster in the project solver.

### 7.3 Mid-size case: complex 2D Laplacian, n = 900

| Phase | Project (ms) | SuperLU (ms) | Slowdown | Project share | Gap (ms) | Gap share |
|---|---:|---:|---:|---:|---:|---:|
| 1. Ordering | 2.870 | 0.174 | 16.5x | 84.2% | 2.696 | 106.0% |
| 2. Symbolic (+ permute) | 0.188 | 0.033 | 5.7x | 5.5% | 0.155 | 6.1% |
| 3. Numeric | 0.298 | 0.635 | 0.47x | 8.7% | -0.337 | -13.2% |
| 4. Triangular solve | 0.052 | 0.036 | 1.4x | 1.5% | 0.016 | 0.6% |
| **Total** | **3.408** | **0.864** | **3.9x** | **100.0%** | **2.544** | **100.0%** |

At this size the project's numeric phase is again faster than SuperLU's, but the ordering gap is already 84% of the runtime.

### 7.4 Scaling of the complex 2D Laplacian

All times are medians of 3 runs. `MD / MMD` uses the direct ordering calls (`get_perm_c` for SuperLU), not the `COLPERM` field inside a full factorization. `Fill increase` is the project's `nnz(L)` relative to SuperLU MMD's Cholesky-equivalent fill.

| `k` | `n` | `nnz(A)` | Project MD (ms) | SuperLU MMD (ms) | MD / MMD | Project fill | MMD fill | Fill increase | Project total (ms) | SuperLU HPD total (ms) | Total slowdown | Project numeric (ms) | SuperLU numeric (ms) | Numeric ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 400 | 1920 | 0.771 | 0.099 | 7.8x | 3729 | 3672 | +1.6% | 0.977 | 0.335 | 2.9x | 0.103 | 0.229 | 0.45x |
| 30 | 900 | 4380 | 2.870 | 0.173 | 16.6x | 10351 | 10098 | +2.5% | 3.408 | 0.864 | 3.9x | 0.298 | 0.635 | 0.47x |
| 40 | 1600 | 7840 | 7.555 | 0.268 | 28.2x | 21504 | 19971 | +7.7% | 8.644 | 1.588 | 5.4x | 0.656 | 1.205 | 0.54x |
| 50 | 2500 | 12300 | 16.317 | 0.403 | 40.5x | 38269 | 35943 | +6.5% | 18.572 | 2.764 | 6.7x | 1.341 | 2.188 | 0.61x |
| 60 | 3600 | 17760 | 31.144 | 0.550 | 56.6x | 60048 | 54215 | +10.8% | 34.245 | 4.203 | 8.1x | 2.163 | 3.462 | 0.62x |
| 80 | 6400 | 31680 | 90.379 | 0.839 | 107.7x | 127600 | 110462 | +15.5% | 98.028 | 8.591 | 11.4x | 5.668 | 7.323 | 0.77x |
| 100 | 10000 | 49600 | 215.973 | 1.298 | 166.4x | 221032 | 185673 | +19.0% | 230.607 | 15.963 | 14.4x | 11.573 | 13.077 | 0.89x |
| 120 | 14400 | 71520 | 401.197 | 1.766 | 227.2x | 342036 | 289179 | +18.3% | 425.479 | 27.416 | 15.5x | 19.114 | 23.816 | 0.80x |

The ordering slowdown grows from 7.8x at `n = 400` to 227x at `n = 14400`. The project's fill is within a few percent of SuperLU MMD at small sizes but becomes 15-19% larger at `n >= 6400`. Despite the larger fill, the project's numeric phase remains faster in absolute wall-clock time because Cholesky does less arithmetic than LU.

### 7.5 Random sparse matrices

| Matrix | `n` | `nnz(A)` | Project MD (ms) | SuperLU MMD (ms) | MD / MMD | Project fill | MMD fill | Project total (ms) | SuperLU HPD total (ms) | Total slowdown | Project numeric (ms) | SuperLU numeric (ms) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `random-60` | 60 | 350 | 0.085 | 0.028 | 3.0x | 438 | 437 | 0.123 | 0.064 | 1.9x | 0.017 | 0.026 |
| `random-300` | 300 | 2198 | 7.091 | 0.626 | 11.3x | 8614 | 8504 | 7.734 | 1.315 | 5.9x | 0.431 | 0.631 |
| `random-900` | 900 | 2628 | 0.963 | 0.220 | 4.4x | 3131 | 3134 | 1.212 | 0.489 | 2.5x | 0.090 | 0.252 |

The random cases show the same qualitative result: MD is 3-11x slower than SuperLU MMD and dominates the total runtime, while numeric factorization is faster in the project solver. The fill quality is essentially identical to MMD on these random matrices.

### 7.6 Ordering quality cross-evaluation

To separate ordering quality from ordering implementation, every ordering was evaluated through the project's own `permute` and `symbolicFactorization`. The number shown is `nnz(L)` from the project's symbolic factorization under that ordering. All three orderings therefore use identical symbolic/factorization code, and the fill difference is purely a property of the ordering.

| `n` | Project MD time (ms) | MMD time (ms) | MD time ratio | Project MD fill | MMD fill | Fill increase | COLAMD time (ms) | COLAMD fill |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 400 | 0.771 | 0.099 | 7.8x | 3729 | 3672 | +1.6% | 0.163 | 5024 |
| 900 | 2.870 | 0.173 | 16.6x | 10351 | 10098 | +2.5% | 0.316 | 15169 |
| 1600 | 7.555 | 0.268 | 28.2x | 21504 | 19971 | +7.7% | 0.470 | 32451 |
| 2500 | 16.317 | 0.403 | 40.5x | 38269 | 35943 | +6.5% | 0.688 | 58252 |
| 3600 | 31.144 | 0.550 | 56.6x | 60048 | 54215 | +10.8% | 0.870 | 89860 |
| 6400 | 90.379 | 0.839 | 107.7x | 127600 | 110462 | +15.5% | 1.268 | 185141 |
| 10000 | 215.973 | 1.298 | 166.4x | 221032 | 185673 | +19.0% | 1.811 | 322875 |
| 14400 | 401.197 | 1.766 | 227.2x | 342036 | 289179 | +18.3% | 2.501 | 503830 |

This table is the key to the diagnosis:

- The project's MD produces fill close to SuperLU MMD at small sizes and 15-19% more fill at large sizes. The ordering quality is therefore not catastrophic, but it degrades with problem size.
- COLAMD, the default SuperLU ordering for general unsymmetric problems, produces much more fill than MMD on these HPD matrices. This is why the stock SuperLU configuration is slower than the HPD no-pivot configuration.
- The MD **time** is 8-227x larger than MMD while the fill is only 1.6-19% larger. Most of the ordering gap is implementation cost, not ordering quality.

### 7.7 Same-ordering numeric comparison

In this experiment the project's MD permutation is fed to SuperLU through `ColPerm = MY_PERMC`, so both solvers use exactly the same permutation. The numeric phase then differs only by the factorization algorithm and kernel.

| `n` | Project numeric (ms) | SuperLU same-order numeric (ms) | Ratio | Project flops (million) | SuperLU flops (million) | Project GFLOP/s | SuperLU GFLOP/s |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 900 | 0.298 | 0.663 | 0.45x | 0.71 | 1.54 | 2.39 | 2.32 |
| 2500 | 1.341 | 2.460 | 0.55x | 4.77 | 9.97 | 3.56 | 4.05 |
| 6400 | 5.668 | 9.343 | 0.61x | 27.64 | 56.72 | 4.88 | 6.07 |
| 10000 | 11.573 | 17.539 | 0.66x | 59.44 | 121.38 | 5.14 | 6.92 |
| 14400 | 19.114 | 29.377 | 0.65x | 103.92 | 211.74 | 5.44 | 7.21 |

On identical fill, the project's numeric phase is 35-55% faster in wall-clock time than SuperLU's numeric phase for these problems. SuperLU has a higher flop rate (about 25-35% higher at `n >= 6400`), but it performs roughly twice as many flops because it computes an LU factorization instead of a Cholesky factorization. The two effects nearly cancel; the project's numeric phase is not the primary bottleneck.

### 7.8 Accuracy

All 12 matrices and all solver modes completed successfully with `info = 0`. The project solver's relative residual is in the range `1.03e-16` to `8.68e-16`; SuperLU HPD is in the range `5.93e-17` to `5.63e-16`. The largest forward errors occur at the largest problem sizes:

| Case | Project rel. residual | Project forward error | SuperLU HPD rel. residual | SuperLU HPD forward error |
|---|---:|---:|---:|---:|
| `laplacian2d-30`, n = 900 | 6.98e-16 | 2.25e-12 | 3.23e-16 | 1.02e-12 |
| `custom-grid-50`, n = 2500 | 5.00e-16 | 3.93e-12 | 4.04e-16 | 3.61e-12 |
| `custom-grid-80`, n = 6400 | 8.68e-16 | 1.66e-11 | 5.63e-16 | 1.20e-11 |
| `custom-grid-100`, n = 10000 | 8.40e-16 | 2.47e-11 | 4.79e-16 | 1.59e-11 |
| `custom-grid-120`, n = 14400 | 7.77e-16 | 3.71e-11 | 5.08e-16 | 2.27e-11 |

The residual is at machine precision in every case. The forward error grows with problem size and conditioning, as expected; both solvers show the same qualitative behavior.

## 8. Bottleneck analysis

### 8.1 What dominates the project runtime

For the complex 2D Laplacian the phase shares at the largest size are:

| Phase | n = 900 | n = 2500 | n = 6400 | n = 10000 | n = 14400 |
|---|---:|---:|---:|---:|---:|
| Ordering | 84.2% | 87.9% | 92.2% | 93.7% | 94.3% |
| Symbolic (+ permute) | 5.5% | 2.8% | 1.3% | 0.9% | 0.7% |
| Numeric | 8.7% | 7.2% | 5.8% | 5.0% | 4.5% |
| Triangular solve | 1.5% | 1.0% | 0.6% | 0.5% | 0.4% |

Ordering is the dominant cost at every size and grows toward 94% of the runtime. Symbolic, numeric, and solve together are less than 6% at `n = 14400`.

### 8.2 What dominates the gap to SuperLU

At `n = 14400`:

| Phase | Gap (ms) | Share of total gap |
|---|---:|---:|
| Ordering | 399.330 | 100.3% |
| Symbolic | 1.955 | 0.5% |
| Numeric | -4.701 | -1.2% |
| Triangular solve | 0.782 | 0.2% |
| **Total** | **398.063** | **100.0%** |

Ordering accounts for essentially all of the total gap. Symbolic and solve contribute less than 1%. Numeric contributes a negative gap because the project's numeric phase is faster than SuperLU's on this problem.

### 8.3 Verdict

**The bottleneck is the implementation of the Minimum Degree ordering, not the numeric factorization.**

The evidence is:

1. **Ordering time:** MD is 3-11x slower than SuperLU MMD on the random matrices and 7.8-227x slower on the 2D Laplacians, reaching 401 ms versus 1.77 ms at `n = 14400`.
2. **Ordering quality:** the fill produced by MD is within 1.6-2.5% of MMD at small sizes and 15-19% larger at large sizes. Even if the fill were reduced to MMD's level, the numeric phase would save at most about 20% of its time, which is less than 4 ms at `n = 14400`.
3. **Numeric kernel:** with the same permutation, the project's numeric phase is 35-55% faster in wall-clock time than SuperLU's LU numeric phase. SuperLU's kernel has a 25-35% higher flop rate at large sizes, but its LU arithmetic is roughly twice as large as Cholesky's, so the kernel efficiency gap does not translate into a total-time problem.
4. **Symbolic and solve:** together they are below 1.5% of the runtime and below 1% of the gap at the largest size.

### 8.4 Root causes in the current MD implementation

The project's MD implementation in `src/md.cpp` is intentionally simple:

- `MD_Graph` stores each vertex adjacency list in `std::set<std::size_t>`. Every insert, erase, and membership test costs `O(log d)` and allocates tree nodes.
- `eliminateVertex` explicitly materializes all fill edges by iterating over every pair of neighbors of the eliminated vertex and inserting both directions into the sets. This is the standard clique update, but it produces a large number of small heap allocations and logarithmic operations.
- `findMinimumDegree` scans all `n` vertices on every elimination step. This alone is `O(n^2)` over the whole ordering, independent of fill.
- The algorithm is a plain minimum-degree method, not a multiple-minimum-degree or approximate-minimum-degree method. It does not use a quotient graph, supernodes, mass elimination, or an external-degree approximation.

SuperLU's MMD implementation uses a quotient-graph representation, multiple elimination, degree lists, and careful array-based data structures. It avoids explicit fill-edge insertion and avoids the linear minimum search. Those are the reasons for the large constant-factor and asymptotic difference.

### 8.5 Recommended optimization order

1. **Replace the MD ordering first.** Use a quotient-graph or AMD-style implementation, or reuse an existing AMD/MMD/COLAMD library. This is the only change that can remove the dominant 94% of runtime. The benchmark's ordering cross-evaluation can be used as a regression test: compare `nnz(L)` and ordering time against SuperLU MMD.
2. **If MD must remain in-house, replace `std::set` with array-based adjacency and degree buckets.** Use linked lists or compressed adjacency with lazy deletion, maintain a bucket queue of degrees, and avoid scanning all vertices for the minimum. This should remove most of the 8-227x constant-factor gap.
3. **Improve ordering quality.** Reducing fill by 15-19% would reduce symbolic and numeric cost proportionally. At `n = 14400`, numeric would fall from about 19 ms to about 16 ms and symbolic from about 3 ms to about 2.5 ms. This is useful but secondary to fixing the ordering time.
4. **Optimize symbolic and triangular solve only after ordering.** Their combined cost is below 1.5% of the total at the largest size. Even eliminating them entirely would change the total time by less than 2%.
5. **Do not prioritize a supernodal BLAS3 numeric kernel for these sizes.** The current scalar Cholesky is already faster in wall-clock time than SuperLU's LU on identical fill. SuperLU's per-flop advantage is real but modest, and LU performs more work. A BLAS3 kernel would be relevant only after the ordering bottleneck is removed and only for much larger or denser factors.

## 9. Comparison with stock SuperLU

The phase-aligned reference is SuperLU's HPD no-pivot mode, because it is the closest analogue of a Cholesky factorization. For completeness, the stock configuration (COLAMD + threshold partial pivoting) was also measured:

| Case | Project total (ms) | SuperLU stock total (ms) | Stock slowdown | SuperLU HPD total (ms) | HPD slowdown |
|---|---:|---:|---:|---:|---:|
| `laplacian2d-30`, n = 900 | 3.408 | 1.280 | 2.7x | 0.864 | 3.9x |
| `custom-grid-50`, n = 2500 | 18.572 | 4.568 | 4.1x | 2.764 | 6.7x |
| `custom-grid-80`, n = 6400 | 98.028 | 14.542 | 6.7x | 8.591 | 11.4x |
| `custom-grid-100`, n = 10000 | 230.607 | 27.565 | 8.4x | 15.963 | 14.4x |
| `custom-grid-120`, n = 14400 | 425.479 | 46.155 | 9.2x | 27.416 | 15.5x |

Stock SuperLU is slower than HPD mode here because COLAMD plus partial pivoting produces much more fill on these HPD matrices. For example, at `random-300` the stock mode has 41620 factor nonzeros versus 16708 in HPD mode; at `n = 14400` it has 1092476 versus 563958. The HPD mode is the fair reference for a Cholesky solver, while the stock result shows what a user gets from SuperLU without tuning for a Hermitian problem.

## 10. Accuracy and validation summary

- The project program was rebuilt with the benchmark changes and still prints `ALL OK`.
- The project solver's own test suite covers the 4x4 tridiagonal matrix, random HPD `n = 60` and `n = 300`, and rejection of non-Hermitian, indefinite, and non-real-diagonal inputs.
- All benchmark matrices are HPD by construction and solve to machine-precision residual.
- SuperLU returned `info = 0` on every measured run. No warning code `n + 1` occurred.
- The relative residual is below `1e-15` for both solvers in every case.
- The largest forward error is `3.71e-11` (project) and `2.27e-11` (SuperLU HPD) at `n = 14400`, consistent with the conditioning of the largest matrix.

## 11. Caveats and limitations

1. **LU versus Cholesky.** SuperLU computes `Pr A Pc = L U`; the project computes `PAP^H = L L^H`. LU generally stores roughly twice as many factor entries and performs roughly twice as many arithmetic operations. Raw `nnz(L)` and raw flop counts are not directly comparable. The report uses the Cholesky-equivalent fill `(nnz(L) + nnz(U)) / 2` and reports effective GFLOP/s separately.
2. **Estimated symbolic split.** SuperLU does not expose a public symbolic-only timer. Phase 2 is estimated as `ETREE + (FACT_first - FACT_reuse)`. The reuse run is warm and preallocated, so the estimate can overstate symbolic time and understate numeric time. The estimate is used only in the no-pivot HPD mode.
3. **Different ordering algorithms.** SuperLU MMD_AT_PLUS_A is not the same algorithm as the project's MD. It is the closest available reference and is evaluated on the same matrix. COLAMD is also reported; it is a general unsymmetric ordering and is not expected to be competitive on these HPD matrices.
4. **Timer resolution at small sizes.** Cases below roughly `n = 1000` have sub-millisecond phases. Scheduler noise, timer resolution, and allocator behavior can dominate. The large 2D cases are the meaningful performance evidence.
5. **Single-threaded BLAS.** OpenBLAS is forced to one thread. Multi-threaded BLAS would change SuperLU's numeric performance and could make its BLAS3 advantage more visible; it would also make the comparison less controlled. The project solver is serial.
6. **Homebrew build.** SuperLU and OpenBLAS are the Homebrew 7.0.1 and 0.3.34 builds. Different BLAS implementations, compilers, or build flags can change absolute times.
7. **Machine specificity.** All numbers are from an Apple M1 with 8 GiB RAM and macOS 27.2. They should not be generalized to other hardware without re-running the benchmark.
8. **No dense HPD validation in the benchmark.** The benchmark uses known HPD constructions and checks the final residual. The `validateHPD` routine is exercised by `./build/program`, not by the timed benchmark, because its dense `O(n^3)` Cholesky validation would distort the measurement.
9. **Warm-up policy.** SuperLU gets one unmeasured warm-up per configuration; the project solver does not. The project solver has no external library initialization, but its first run may still benefit from allocator and cache warmth. The median of three runs is used to reduce this effect.
10. **Phase 4 definition.** The triangular solve includes the permutation of the right-hand side and the final permutation of the solution, matching the paper's phase 4.
11. **Generated matrix suite.** The benchmark uses the project's generated HPD matrices and structured 2D Laplacians, not external SuiteSparse matrices. Real application matrices can have different degree distributions and fill behavior, so the ordering bottleneck should be re-checked on representative application data.

## 12. All solver runs

The table below contains every solver run in the final benchmark, including the stock SuperLU mode and the same-order numeric diagnostic. Times are milliseconds (median of 3). `fill` is `nnz(factor) / nnz(A)`; for SuperLU it uses `nnz(L) + nnz(U) - n`, while for the project it uses `nnz(L)`. The same-order rows have zero ordering time because the project's permutation is supplied. Forward errors are summarized in section 7.8; the simplified CSV written by `--csv` contains one row per matrix and phase with the ours-minus-SuperLU time difference in milliseconds and the ours-over-SuperLU ratio.

| Matrix | n | nnz(A) | Solver/mode | Order | Symbolic | Numeric | Solve | Total | nnz factor | fill | rel. residual |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| tridiagonal-300 | 300 | 898 | ours/md | 0.066 | 0.035 | 0.016 | 0.008 | 0.125 | 599 | 0.67 | 1.03e-16 |
| tridiagonal-300 | 300 | 898 | superlu/default | 0.022 | 0.008 | 0.032 | 0.008 | 0.070 | 990 | 1.10 | 5.97e-17 |
| tridiagonal-300 | 300 | 898 | superlu/hpd | 0.015 | 0.005 | 0.033 | 0.005 | 0.058 | 898 | 1.00 | 5.93e-17 |
| tridiagonal-300 | 300 | 898 | superlu/same-order | 0.000 | 0.005 | 0.035 | 0.008 | 0.048 | 988 | 1.10 | 5.97e-17 |
| random-60 | 60 | 350 | ours/md | 0.085 | 0.019 | 0.017 | 0.003 | 0.123 | 438 | 1.25 | 1.17e-16 |
| random-60 | 60 | 350 | superlu/default | 0.027 | 0.004 | 0.038 | 0.004 | 0.072 | 1662 | 4.75 | 1.26e-16 |
| random-60 | 60 | 350 | superlu/hpd | 0.029 | 0.006 | 0.026 | 0.002 | 0.064 | 828 | 2.37 | 7.81e-17 |
| random-60 | 60 | 350 | superlu/same-order | 0.000 | 0.003 | 0.021 | 0.002 | 0.026 | 816 | 2.33 | 1.06e-16 |
| random-300 | 300 | 2198 | ours/md | 7.091 | 0.179 | 0.431 | 0.042 | 7.734 | 8614 | 3.92 | 3.54e-16 |
| random-300 | 300 | 2198 | superlu/default | 0.339 | 0.073 | 2.243 | 0.050 | 2.646 | 41620 | 18.94 | 2.65e-16 |
| random-300 | 300 | 2198 | superlu/hpd | 0.623 | 0.037 | 0.631 | 0.022 | 1.315 | 16708 | 7.60 | 1.86e-16 |
| random-300 | 300 | 2198 | superlu/same-order | 0.000 | 0.028 | 0.651 | 0.023 | 0.691 | 16928 | 7.70 | 1.74e-16 |
| laplacian2d-20 | 400 | 1920 | ours/md | 0.771 | 0.085 | 0.103 | 0.018 | 0.977 | 3729 | 1.94 | 2.40e-16 |
| laplacian2d-20 | 400 | 1920 | superlu/default | 0.147 | 0.023 | 0.268 | 0.025 | 0.471 | 11742 | 6.12 | 2.58e-16 |
| laplacian2d-20 | 400 | 1920 | superlu/hpd | 0.081 | 0.013 | 0.229 | 0.016 | 0.335 | 6944 | 3.62 | 3.49e-16 |
| laplacian2d-20 | 400 | 1920 | superlu/same-order | 0.000 | 0.008 | 0.233 | 0.015 | 0.256 | 7058 | 3.68 | 2.13e-16 |
| laplacian2d-30 | 900 | 4380 | ours/md | 2.870 | 0.188 | 0.298 | 0.052 | 3.408 | 10351 | 2.36 | 6.98e-16 |
| laplacian2d-30 | 900 | 4380 | superlu/default | 0.308 | 0.043 | 0.868 | 0.064 | 1.280 | 35150 | 8.03 | 2.94e-16 |
| laplacian2d-30 | 900 | 4380 | superlu/hpd | 0.174 | 0.033 | 0.635 | 0.036 | 0.864 | 19296 | 4.41 | 3.23e-16 |
| laplacian2d-30 | 900 | 4380 | superlu/same-order | 0.000 | 0.026 | 0.663 | 0.038 | 0.713 | 19802 | 4.52 | 3.61e-16 |
| random-900 | 900 | 2628 | ours/md | 0.963 | 0.141 | 0.090 | 0.015 | 1.212 | 3131 | 1.19 | 1.26e-16 |
| random-900 | 900 | 2628 | superlu/default | 0.264 | 0.032 | 0.348 | 0.039 | 0.683 | 12968 | 4.93 | 1.57e-16 |
| random-900 | 900 | 2628 | superlu/hpd | 0.197 | 0.028 | 0.252 | 0.017 | 0.489 | 5368 | 2.04 | 1.37e-16 |
| random-900 | 900 | 2628 | superlu/same-order | 0.000 | 0.024 | 0.236 | 0.020 | 0.280 | 5386 | 2.05 | 1.62e-16 |
| laplacian2d-40 | 1600 | 7840 | ours/md | 7.555 | 0.332 | 0.656 | 0.105 | 8.644 | 21504 | 2.74 | 4.88e-16 |
| laplacian2d-40 | 1600 | 7840 | superlu/default | 0.471 | 0.078 | 1.809 | 0.117 | 2.471 | 73600 | 9.39 | 3.36e-16 |
| laplacian2d-40 | 1600 | 7840 | superlu/hpd | 0.254 | 0.051 | 1.205 | 0.070 | 1.588 | 38342 | 4.89 | 3.86e-16 |
| laplacian2d-40 | 1600 | 7840 | superlu/same-order | 0.000 | 0.022 | 1.318 | 0.064 | 1.410 | 41408 | 5.28 | 4.02e-16 |
| custom-grid-50 | 2500 | 12300 | ours/md | 16.317 | 0.513 | 1.341 | 0.194 | 18.572 | 38269 | 3.11 | 5.00e-16 |
| custom-grid-50 | 2500 | 12300 | superlu/default | 0.709 | 0.263 | 3.390 | 0.199 | 4.568 | 130686 | 10.62 | 3.68e-16 |
| custom-grid-50 | 2500 | 12300 | superlu/hpd | 0.388 | 0.055 | 2.188 | 0.111 | 2.764 | 69386 | 5.64 | 4.04e-16 |
| custom-grid-50 | 2500 | 12300 | superlu/same-order | 0.000 | 0.078 | 2.460 | 0.112 | 2.632 | 74038 | 6.02 | 4.67e-16 |
| custom-grid-60 | 3600 | 17760 | ours/md | 31.144 | 0.718 | 2.163 | 0.286 | 34.245 | 60048 | 3.38 | 6.38e-16 |
| custom-grid-60 | 3600 | 17760 | superlu/default | 0.870 | 0.248 | 5.361 | 0.293 | 6.779 | 200258 | 11.28 | 3.63e-16 |
| custom-grid-60 | 3600 | 17760 | superlu/hpd | 0.525 | 0.040 | 3.462 | 0.160 | 4.203 | 104830 | 5.90 | 4.08e-16 |
| custom-grid-60 | 3600 | 17760 | superlu/same-order | 0.000 | 0.059 | 3.834 | 0.159 | 4.051 | 116496 | 6.56 | 3.76e-16 |
| custom-grid-80 | 6400 | 31680 | ours/md | 90.379 | 1.310 | 5.668 | 0.592 | 98.028 | 127600 | 4.03 | 8.68e-16 |
| custom-grid-80 | 6400 | 31680 | superlu/default | 1.301 | 0.181 | 12.399 | 0.572 | 14.542 | 407446 | 12.86 | 4.63e-16 |
| custom-grid-80 | 6400 | 31680 | superlu/hpd | 0.855 | 0.119 | 7.323 | 0.298 | 8.591 | 214524 | 6.77 | 5.63e-16 |
| custom-grid-80 | 6400 | 31680 | superlu/same-order | 0.000 | 0.152 | 9.343 | 0.325 | 9.838 | 248800 | 7.85 | 3.83e-16 |
| custom-grid-100 | 10000 | 49600 | ours/md | 215.973 | 1.976 | 11.573 | 1.085 | 230.607 | 221032 | 4.46 | 8.40e-16 |
| custom-grid-100 | 10000 | 49600 | superlu/default | 1.793 | 0.249 | 24.480 | 1.032 | 27.565 | 704340 | 14.20 | 4.40e-16 |
| custom-grid-100 | 10000 | 49600 | superlu/hpd | 1.285 | 0.453 | 13.077 | 0.494 | 15.963 | 361346 | 7.29 | 4.79e-16 |
| custom-grid-100 | 10000 | 49600 | superlu/same-order | 0.000 | 0.140 | 17.539 | 0.545 | 18.206 | 432064 | 8.71 | 4.85e-16 |
| custom-grid-120 | 14400 | 71520 | ours/md | 401.197 | 2.954 | 19.114 | 1.626 | 425.479 | 342036 | 4.78 | 7.77e-16 |
| custom-grid-120 | 14400 | 71520 | superlu/default | 2.409 | 0.734 | 41.413 | 1.650 | 46.155 | 1092476 | 15.28 | 4.80e-16 |
| custom-grid-120 | 14400 | 71520 | superlu/hpd | 1.867 | 0.998 | 23.816 | 0.844 | 27.416 | 563958 | 7.89 | 5.08e-16 |
| custom-grid-120 | 14400 | 71520 | superlu/same-order | 0.000 | 0.503 | 29.377 | 0.867 | 30.693 | 669672 | 9.36 | 6.00e-16 |

## 13. Files added for the benchmark

| File | Role |
|---|---|
| `include/generators.hpp` | Shared random HPD, tridiagonal HPD, and complex 2D Laplacian generators |
| `benchmark/superlu_adapter.hpp` | SuperLU-facing result and option declarations |
| `benchmark/superlu_adapter.cpp` | CSC conversion, `get_perm_c`, `zgssvx`, phase timing, same-pattern refactorization |
| `benchmark/bench_superlu.cpp` | Matrix suite, four-phase timing, accuracy checks, ordering cross-evaluation, bottleneck diagnosis |
| `benchmark/run.sh` | Single-thread BLAS environment and benchmark launcher |
| `run-benchmark.sh` | One-command build and full benchmark run on macOS/Linux without downloading dependencies |
| `run-benchmark.ps1` | Native Windows PowerShell build and benchmark launcher |
| `run-benchmark.bat` | Windows Command Prompt wrapper for `run-benchmark.ps1` |
| `CMakeLists.txt` | Cross-platform CMake build used by Windows and optional elsewhere |
| `Makefile` | `program` and `bench` targets, Homebrew/pkg-config/Linux prefix detection |
| `docs/running.md` | Per-platform benchmark instructions |
| `docs/benchmark.md` | This document |

The project solver itself is unchanged in its algorithmic behavior. The only source change outside the benchmark is the removal of the duplicated `randomHPD` generator from `src/main.cpp`; it now uses the shared implementation in `include/generators.hpp`. The solver still passes its original test suite with `ALL OK`.
