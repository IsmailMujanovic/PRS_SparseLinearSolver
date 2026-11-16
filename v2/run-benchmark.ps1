[CmdletBinding()]
param(
    [string]$SuperLURoot = $env:SUPERLU_ROOT,
    [string]$BuildDir = "build-windows",
    [string]$CMakePath = "",
    [int]$Repeat = 0,
    [switch]$Large,
    [switch]$Detailed,
    [int[]]$Grid = @(),
    [string]$Csv = "",
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$BenchmarkArgs
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $Root

function Find-CMake {
    if ($CMakePath) {
        if (Test-Path $CMakePath) {
            return (Resolve-Path $CMakePath).Path
        }
        throw "CMake was not found at $CMakePath"
    }

    $command = Get-Command cmake -ErrorAction SilentlyContinue
    if ($command) {
        return $command.Source
    }

    $candidates = @()
    if ($env:ProgramFiles) {
        $candidates += (Join-Path $env:ProgramFiles "CMake\bin\cmake.exe")
    }
    if (${env:ProgramFiles(x86)}) {
        $candidates += (Join-Path ${env:ProgramFiles(x86)} "CMake\bin\cmake.exe")
    }
    foreach ($candidate in $candidates) {
        if ($candidate -and (Test-Path $candidate)) {
            return $candidate
        }
    }

    throw "CMake was not found. Install CMake or pass -CMakePath."
}

$cmake = Find-CMake

if (-not $SuperLURoot -and $env:VCPKG_ROOT) {
    $candidate = Join-Path $env:VCPKG_ROOT "installed\x64-windows"
    if (Test-Path $candidate) {
        $SuperLURoot = $candidate
    }
}

if ($SuperLURoot) {
    if (-not (Test-Path $SuperLURoot)) {
        throw "SuperLU root was not found: $SuperLURoot"
    }
    $SuperLURoot = (Resolve-Path $SuperLURoot).Path
    $headerA = Join-Path $SuperLURoot "include\slu_zdefs.h"
    $headerB = Join-Path $SuperLURoot "include\superlu\slu_zdefs.h"
    if (-not (Test-Path $headerA) -and -not (Test-Path $headerB)) {
        Write-Warning "slu_zdefs.h was not found under $SuperLURoot\include. CMake may still find SuperLU through its package configuration."
    }
}

if ([System.IO.Path]::IsPathRooted($BuildDir)) {
    $BuildPath = $BuildDir
} else {
    $BuildPath = Join-Path $Root $BuildDir
}

$configureArgs = @("-S", $Root, "-B", $BuildPath, "-DCMAKE_BUILD_TYPE=Release")
if ($SuperLURoot) {
    $configureArgs += "-DSUPERLU_ROOT=$SuperLURoot"
}

$vcpkgRoot = $null
if ($SuperLURoot -match '[\\/]installed[\\/]') {
    $parts = $SuperLURoot -split '[\\/]'
    $installedIndex = [Array]::LastIndexOf($parts, "installed")
    if ($installedIndex -gt 0) {
        $vcpkgRoot = ($parts[0..($installedIndex - 1)] -join "\")
    }
} elseif (-not $SuperLURoot -and $env:VCPKG_ROOT) {
    $vcpkgRoot = $env:VCPKG_ROOT
}

if ($vcpkgRoot) {
    $toolchain = Join-Path $vcpkgRoot "scripts\buildsystems\vcpkg.cmake"
    if (Test-Path $toolchain) {
        $configureArgs += "-DCMAKE_TOOLCHAIN_FILE=$toolchain"
        if ($SuperLURoot) {
            $triplet = Split-Path $SuperLURoot -Leaf
            if ($triplet -match '^(x64|x86|arm64)-windows') {
                $configureArgs += "-DVCPKG_TARGET_TRIPLET=$triplet"
            }
        }
    }
}

& $cmake @configureArgs
if ($LASTEXITCODE -ne 0) {
    throw "CMake configuration failed. Install SuperLU and pass -SuperLURoot, or set SUPERLU_ROOT."
}

& $cmake --build $BuildPath --config Release --target bench
if ($LASTEXITCODE -ne 0) {
    throw "Building the benchmark failed."
}

$executable = Get-ChildItem -Path $BuildPath -Filter "bench.exe" -Recurse -File |
    Sort-Object LastWriteTime -Descending |
    Select-Object -First 1
if (-not $executable) {
    throw "bench.exe was not produced under $BuildPath."
}

$env:OPENBLAS_NUM_THREADS = "1"
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"

$arguments = @()
if ($BenchmarkArgs) {
    $arguments += $BenchmarkArgs
}
if ($Repeat -gt 0) {
    $arguments += @("--repeat", "$Repeat")
}
if ($Large) {
    $arguments += "--large"
}
if ($Detailed) {
    $arguments += "--detailed"
}
foreach ($gridSize in $Grid) {
    $arguments += @("--grid", "$gridSize")
}
if ($Csv) {
    $arguments += @("--csv", $Csv)
}

$csvPath = $null
if ($arguments.Count -eq 0) {
    $arguments = @(
        "--repeat", "3",
        "--large",
        "--grid", "50",
        "--grid", "60",
        "--grid", "80",
        "--grid", "100",
        "--grid", "120"
    )
} elseif ($Csv) {
    $csvPath = $Csv
}

& $executable.FullName @arguments
$exitCode = $LASTEXITCODE

if ($csvPath) {
    Write-Host ""
    Write-Host "Results CSV: $csvPath"
}

exit $exitCode
