# Running the benchmark on each machine

The benchmark can be run on macOS, Linux, and native Windows. The dependency on SuperLU and a BLAS library is the same everywhere; only the installation and launcher differ.

The launcher scripts never download anything. They build the benchmark from the source in this repository and run it. SuperLU and the BLAS library must already be installed on the machine.

## Quick reference

| Platform | Launcher | One-time prerequisites |
|---|---|---|
| macOS | `./run-benchmark.sh` | Xcode Command Line Tools, Homebrew, `superlu` |
| Linux | `./run-benchmark.sh` | C++ compiler, make, SuperLU development package, BLAS |
| Windows | `run-benchmark.bat` | Visual Studio C++ Build Tools, CMake, SuperLU built for MSVC |
| Windows PowerShell | `.\run-benchmark.ps1` | Same as above |
| WSL | `./run-benchmark.sh` | Same as Linux |

All commands are run from the project root.

## Common behavior

The launchers:

- build the benchmark if necessary;
- run the full benchmark suite by default (`--repeat 3`, 2D Laplacians up to `n = 14400`);
- set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, and `MKL_NUM_THREADS=1`;
- print a simple four-phase comparison with only the project time, the SuperLU time, and the total;
- write a CSV only when `--csv FILE` (Windows: `-Csv FILE`) is explicitly requested.

The benchmark can also be run directly after building:

```bash
./build/bench --repeat 3
```

On Windows:

```powershell
.\build-windows\Release\bench.exe --repeat 3
```

The complete methodology, results, and analysis are in [benchmark.md](benchmark.md).

## macOS

### Prerequisites

Install Xcode Command Line Tools and Homebrew, then install SuperLU:

```bash
xcode-select --install
brew install superlu
```

Homebrew installs OpenBLAS automatically as a SuperLU dependency.

### Run

```bash
./run-benchmark.sh
```

Quick run:

```bash
./run-benchmark.sh --repeat 1
```

Custom run:

```bash
./run-benchmark.sh --repeat 3 --large --grid 50 --grid 100
```

### Manual build

```bash
make program
make bench
./build/program
./build/bench --repeat 3
```

The Makefile finds Homebrew SuperLU and OpenBLAS through `brew --prefix`.

### CMake build

CMake is optional on macOS; the Makefile is sufficient. To use CMake:

```bash
cmake -S . -B build-cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DSUPERLU_ROOT="$(brew --prefix superlu)"
cmake --build build-cmake --config Release --target bench
./build-cmake/bench --repeat 3
```

## Linux

### Debian / Ubuntu

```bash
sudo apt update
sudo apt install build-essential pkg-config libsuperlu-dev libopenblas-dev
```

If you want to use the CMake build as well:

```bash
sudo apt install cmake
```

### Fedora

```bash
sudo dnf install gcc-c++ make pkgconf-pkg-config superlu-devel openblas-devel
```

### Arch Linux

```bash
sudo pacman -S base-devel superlu openblas
```

### Run

```bash
./run-benchmark.sh
```

The Makefile uses `pkg-config` when available and also recognizes the Debian header layout `/usr/include/superlu/`. If the automatic detection picks the wrong installation, override the paths:

```bash
make bench \
  SUPERLU_CFLAGS="-I/usr/include/superlu" \
  SUPERLU_LIBS="-lsuperlu" \
  BLAS_LIBS="-lopenblas"
```

Then run:

```bash
./benchmark/run.sh --repeat 3
```

### CMake build on Linux

```bash
cmake -S . -B build-cmake -DCMAKE_BUILD_TYPE=Release -DSUPERLU_ROOT=/usr
cmake --build build-cmake --config Release --target bench
./build-cmake/bench --repeat 3
```

If CMake cannot find the package configuration directly, set its directory explicitly:

```bash
cmake -S . -B build-cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -Dsuperlu_DIR=/usr/lib/x86_64-linux-gnu/cmake/superlu
```

## Windows

Native Windows support uses CMake and PowerShell. The benchmark source is portable C++17; the only platform-specific parts are the launcher and the SuperLU library build, which must use the same compiler toolchain as the benchmark.

### Prerequisites

1. **Visual Studio Build Tools 2022** (or Visual Studio 2022) with the **Desktop development with C++** workload.
2. **CMake**, either on `PATH` or in its default install location.
3. **SuperLU 7.0.1 built for MSVC.**

The easiest way to get a matching SuperLU build is vcpkg. The vcpkg setup is a one-time download; the benchmark launcher itself does not download anything:

```bat
git clone https://github.com/microsoft/vcpkg C:\vcpkg
C:\vcpkg\bootstrap-vcpkg.bat
C:\vcpkg\vcpkg install superlu:x64-windows
```

The installed tree is:

```text
C:\vcpkg\installed\x64-windows
```

Set it for the current PowerShell session:

```powershell
$env:SUPERLU_ROOT = "C:\vcpkg\installed\x64-windows"
```

Alternatively, pass it explicitly with `-SuperLURoot`. When the SuperLU path is inside a vcpkg installation, or when `VCPKG_ROOT` is set, the launcher also passes the vcpkg CMake toolchain file so that SuperLU's BLAS dependency resolves correctly.

### Run from PowerShell

```powershell
.\run-benchmark.ps1 -SuperLURoot C:\vcpkg\installed\x64-windows
```

If `SUPERLU_ROOT` or `VCPKG_ROOT` is set, the path argument is optional:

```powershell
.\run-benchmark.ps1
```

Quick run:

```powershell
.\run-benchmark.ps1 -SuperLURoot C:\vcpkg\installed\x64-windows -Repeat 1
```

### Run from Command Prompt

```bat
run-benchmark.bat -SuperLURoot C:\vcpkg\installed\x64-windows
```

The `.bat` file bypasses the PowerShell execution policy for the launcher. You can also run the PowerShell file directly:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\run-benchmark.ps1 -SuperLURoot C:\vcpkg\installed\x64-windows
```

### CMake build manually

```powershell
cmake -S . -B build-windows -DSUPERLU_ROOT=C:\vcpkg\installed\x64-windows
cmake --build build-windows --config Release --target bench
.\build-windows\Release\bench.exe --repeat 3
```

If CMake reports that no compiler was found, run the commands from **x64 Native Tools Command Prompt for VS 2022** or from **Developer PowerShell for VS 2022**.

If CMake finds SuperLU headers and the library but fails to link BLAS symbols, pass the BLAS library explicitly:

```powershell
cmake -S . -B build-windows `
  -DSUPERLU_ROOT=C:\vcpkg\installed\x64-windows `
  -DSUPERLU_BLAS_LIBRARIES=C:\vcpkg\installed\x64-windows\lib\openblas.lib
```

### WSL

If native Windows setup is inconvenient, install WSL and use the Linux instructions inside it:

```bash
sudo apt install build-essential pkg-config libsuperlu-dev libopenblas-dev
./run-benchmark.sh
```

WSL results are Linux results, not native Windows results.

## Options

The macOS/Linux launcher forwards the benchmark options directly. The Windows launcher exposes the common options as PowerShell parameters and forwards any remaining arguments to the benchmark.

| macOS/Linux | Windows PowerShell | Meaning |
|---|---|---|
| `--repeat N` | `-Repeat N` | Measured repetitions; the median is reported. Default: 3. |
| `--large` | `-Large` | Add `random-900` and `laplacian2d-40`. |
| `--grid K` | `-Grid K` or `-Grid K1,K2` | Add one or more `K x K` complex 2D Laplacians. |
| `--random N:DENSITY:SEED` | `-BenchmarkArgs @("--random","N:DENSITY:SEED")` | Add a random HPD matrix. |
| `--detailed` | `-Detailed` | Print the extra diagnostic tables and accuracy information. |
| `--csv FILE` | `-Csv FILE` | Write every solver run to a CSV file. |
| `--help` | `-BenchmarkArgs @("--help")` | Print benchmark usage. |

Examples:

```bash
./run-benchmark.sh --repeat 1
./run-benchmark.sh --grid 40 --grid 80
./run-benchmark.sh --random 2000:0.001:42
./run-benchmark.sh --detailed
```

```powershell
.\run-benchmark.ps1 -Repeat 1
.\run-benchmark.ps1 -Grid 40,80
.\run-benchmark.ps1 -Repeat 3 -Large -Grid 50,60,80,100,120
.\run-benchmark.ps1 -Detailed
```

## Output

The default console output is a compact table per matrix:

```text
matrix-name  n=...  nnz(A)=...
Phase                        ours (ms)    SuperLU (ms)   ours - SuperLU (ms)  ours / SuperLU (x)
1 Ordering                       ...             ...                 ...                 ...
2 Symbolic factorization         ...             ...                 ...                 ...
3 Numeric factorization          ...             ...                 ...                 ...
4 Triangular solution            ...             ...                 ...                 ...
Total                            ...             ...                 ...                 ...
```

Only the four phases, the total, the `ours - SuperLU (ms)` difference, and the `ours / SuperLU (x)` ratio are shown. A positive difference means the project phase is slower; a negative difference means it is faster. A ratio above 1 means SuperLU is faster; a ratio below 1 means the project phase is faster. The extra diagnostic tables (SuperLU variants, ordering quality, fill, accuracy, bottleneck diagnosis) are available with `--detailed` or `-Detailed`.

The CSV is written only when explicitly requested:

```bash
./run-benchmark.sh --csv build/benchmark_results.csv
```

```powershell
.\run-benchmark.ps1 -Csv build-windows\benchmark_results.csv
```

The CSV contains one row per matrix and phase:

```text
matrix,n,phase,ours_ms,superlu_ms,ours_minus_superlu_ms,ours_over_superlu_x
```

The difference is in milliseconds: positive means the project phase is slower, negative means it is faster. The final column is the ratio: above 1 means SuperLU is faster, below 1 means the project phase is faster.

## Troubleshooting

### SuperLU not found on macOS or Linux

Check that the development files are installed:

```bash
find /usr /usr/local /opt/homebrew -name slu_zdefs.h 2>/dev/null
```

Then pass the include and library paths explicitly:

```bash
make bench \
  SUPERLU_CFLAGS="-I/path/to/include" \
  SUPERLU_LIBS="-L/path/to/lib -lsuperlu" \
  BLAS_LIBS="-lopenblas"
```

### SuperLU not found on Windows

Verify that the directory contains:

```text
include\slu_zdefs.h
lib\superlu.lib
```

or:

```text
include\superlu\slu_zdefs.h
lib\superlu.lib
```

Then pass that directory with `-SuperLURoot`. If using vcpkg, the default is `C:\vcpkg\installed\x64-windows`.

### CMake cannot find SuperLU on Windows

Pass the package configuration directory explicitly:

```powershell
cmake -S . -B build-windows `
  -Dsuperlu_DIR=C:\vcpkg\installed\x64-windows\share\superlu
```

### Wrong compiler or ABI mismatch

MSVC cannot link a SuperLU library built with MinGW, and vice versa. Build SuperLU and the benchmark with the same toolchain. With vcpkg, use the `x64-windows` triplet for MSVC.

### PowerShell refuses to run the script

Use the `.bat` launcher, or start PowerShell with:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\run-benchmark.ps1
```

### CMake generator errors

Run from a Visual Studio developer prompt. You can also select a generator explicitly, for example:

```powershell
cmake -S . -B build-windows -G "Visual Studio 17 2022" -A x64 -DSUPERLU_ROOT=C:\vcpkg\installed\x64-windows
```

### BLAS threading

The launchers set the BLAS thread environment variables to `1`. If you run the benchmark executable directly, set them yourself for a fair comparison:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

```powershell
$env:OPENBLAS_NUM_THREADS = "1"
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
```

### Timings look different from `docs/benchmark.md`

Expected. The recorded results come from an Apple M1 with Apple clang 21.0.0, SuperLU 7.0.1, and OpenBLAS 0.3.34. Different CPUs, compilers, BLAS builds, and operating systems produce different absolute times. The phase comparison and bottleneck method remain valid, but the numbers should be regenerated on each machine.

## Notes on reproducibility

- The launchers do not download dependencies.
- The benchmark uses only the C++ standard library, SuperLU, and the BLAS implementation that SuperLU uses.
- SuperLU 7.0.1 is the tested version. Other versions may work, but the API has not been verified across every release.
- On Windows, SuperLU's internal phase timer uses `clock()` rather than wall-clock time. The benchmark runs single-threaded, so the difference is small, but Windows phase timings are not perfectly identical to macOS/Linux.
- The full benchmark is intended to run on an otherwise idle machine.
