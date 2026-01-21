@echo off
setlocal

if not exist eagar_tsai_integrand.c (
  echo ERROR: eagar_tsai_integrand.c not found in current directory.
  exit /b 1
)

where cl >NUL 2>&1
if %ERRORLEVEL%==0 (
  echo Building with MSVC...
  cl /O2 /LD eagar_tsai_integrand.c /Fe:libeagar_tsai_integrand.dll
  if %ERRORLEVEL%==0 (
    echo Built libeagar_tsai_integrand.dll with MSVC.
    exit /b 0
  ) else (
    echo MSVC build failed.
  )
)

where gcc >NUL 2>&1
if %ERRORLEVEL%==0 (
  echo Building with MinGW gcc...
  gcc -O3 -shared -o libeagar_tsai_integrand.dll eagar_tsai_integrand.c -Wl,--out-implib,libeagar_tsai_integrand.a
  if %ERRORLEVEL%==0 (
    echo Built libeagar_tsai_integrand.dll with MinGW gcc.
    exit /b 0
  ) else (
    echo MinGW gcc build failed.
  )
)

echo ERROR: No compiler found. Install Visual Studio Build Tools (MSVC) or MinGW-w64.
exit /b 1
