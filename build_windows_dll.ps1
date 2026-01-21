$ErrorActionPreference = "Stop"

if (!(Test-Path "eagar_tsai_integrand.c")) {
  Write-Error "eagar_tsai_integrand.c not found in current directory."
}

$cl = Get-Command cl -ErrorAction SilentlyContinue
if ($cl) {
  Write-Host "Building with MSVC..."
  & cl /O2 /LD eagar_tsai_integrand.c /Fe:libeagar_tsai_integrand.dll
  if ($LASTEXITCODE -eq 0) {
    Write-Host "Built libeagar_tsai_integrand.dll with MSVC."
    exit 0
  } else {
    Write-Warning "MSVC build failed."
  }
}

$gcc = Get-Command gcc -ErrorAction SilentlyContinue
if ($gcc) {
  Write-Host "Building with MinGW gcc..."
  & gcc -O3 -shared -o libeagar_tsai_integrand.dll eagar_tsai_integrand.c -Wl,--out-implib,libeagar_tsai_integrand.a
  if ($LASTEXITCODE -eq 0) {
    Write-Host "Built libeagar_tsai_integrand.dll with MinGW gcc."
    exit 0
  } else {
    Write-Warning "MinGW gcc build failed."
  }
}

Write-Error "No compiler found. Install Visual Studio Build Tools (MSVC) or MinGW-w64."
