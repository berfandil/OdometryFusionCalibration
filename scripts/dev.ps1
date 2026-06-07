# dev.ps1 - one-command build/test gate for the autonomous workflow.
# Encapsulates the VS2022 dev-environment build loop so sub-agents needn't
# re-derive it. Auto-discovers Visual Studio via vswhere.
#
# Usage:
#   powershell -File scripts/dev.ps1 -Task test       # configure (if needed) + build + ctest  [default]
#   powershell -File scripts/dev.ps1 -Task build
#   powershell -File scripts/dev.ps1 -Task configure
#   powershell -File scripts/dev.ps1 -Task clean
param(
    [ValidateSet('configure', 'build', 'test', 'clean')]
    [string]$Task = 'test'
)
$ErrorActionPreference = 'Stop'

$root  = Split-Path $PSScriptRoot -Parent
$build = Join-Path $root 'build'

$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
if (-not (Test-Path $vswhere)) { throw "vswhere not found - is Visual Studio installed?" }
$vs = & $vswhere -latest -property installationPath
if (-not $vs) { throw "No Visual Studio installation found via vswhere." }

$vcvars = Join-Path $vs 'VC\Auxiliary\Build\vcvars64.bat'
$cmExt  = Join-Path $vs 'Common7\IDE\CommonExtensions\Microsoft\CMake'
$cmake  = Join-Path $cmExt 'CMake\bin\cmake.exe'
$ctest  = Join-Path $cmExt 'CMake\bin\ctest.exe'
$ninja  = Join-Path $cmExt 'Ninja\ninja.exe'
foreach ($t in @($vcvars, $cmake, $ninja)) {
    if (-not (Test-Path $t)) { throw "Missing build tool: $t" }
}

# Run a command inside the x64 dev environment.
function Invoke-Dev([string]$inner) {
    cmd /c "`"$vcvars`" >NUL && $inner"
    if ($LASTEXITCODE -ne 0) { throw "Command failed (exit $LASTEXITCODE): $inner" }
}

function Invoke-Configure {
    # OFC_BUILD_ADAPTERS=ON (Slice 13): build the relaxed-edge adapters + their doctest
    # target so the gate covers them alongside the core 'unit' test. NOTE: configure is
    # CACHED below (only re-run when build.ninja is absent), so after changing this flag
    # you must run -Task clean first, then -Task test, for it to take effect.
    Invoke-Dev "`"$cmake`" -G Ninja -S `"$root`" -B `"$build`" -DCMAKE_MAKE_PROGRAM=`"$ninja`" -DCMAKE_BUILD_TYPE=Release -DOFC_BUILD_ADAPTERS=ON"
}

switch ($Task) {
    'clean' {
        if (Test-Path $build) { Remove-Item -Recurse -Force $build }
        Write-Host "Cleaned $build"
    }
    'configure' { Invoke-Configure }
    'build' {
        if (-not (Test-Path (Join-Path $build 'build.ninja'))) { Invoke-Configure }
        Invoke-Dev "`"$cmake`" --build `"$build`""
    }
    'test' {
        if (-not (Test-Path (Join-Path $build 'build.ninja'))) { Invoke-Configure }
        Invoke-Dev "`"$cmake`" --build `"$build`""
        Invoke-Dev "`"$ctest`" --test-dir `"$build`" --output-on-failure"
    }
}
Write-Host "dev.ps1 -Task $Task : OK"
