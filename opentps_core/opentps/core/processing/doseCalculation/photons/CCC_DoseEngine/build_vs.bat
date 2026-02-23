@echo off
setlocal enabledelayedexpansion
REM Build CCC_DoseEngine with Visual Studio (cl.exe).
REM Run from: "x64 Native Tools Command Prompt for VS 2022" (or your VS version),
REM    or run this from a normal cmd - it will try to find and call vcvars64.bat.

set "SCRIPT_DIR=%~dp0"
cd /d "%SCRIPT_DIR%"

REM If cl.exe is not in PATH, try to set up VS environment
where cl.exe >nul 2>&1
if errorlevel 1 (
  echo cl.exe not in PATH. Trying to find Visual Studio...
  for %%V in (2022 2019 2017 18) do (
    for %%E in (Community Professional Enterprise BuildTools) do (
      set "VCVARS=C:\Program Files\Microsoft Visual Studio\%%V\%%E\VC\Auxiliary\Build\vcvars64.bat"
      if exist "!VCVARS!" (
        call "!VCVARS!" >nul 2>&1
        goto :env_ready
      )
    )
  )
  echo ERROR: Visual Studio not found. Please run this from "x64 Native Tools Command Prompt".
  exit /b 1
)
:env_ready

echo Compiling CCC_DoseEngine...
cl /nologo /O2 /EHsc /c convolution.cpp calc_deff.cpp parse_func.cpp make_poly.cpp raytrace.cpp terma_dose_masks.cpp terma_kerma.cpp calc_dose.cpp util.cpp
if errorlevel 1 ( echo Build failed at compile. & exit /b 1 )

link /nologo /OUT:CCC_DoseEngine_win.exe convolution.obj calc_deff.obj parse_func.obj make_poly.obj raytrace.obj terma_dose_masks.obj terma_kerma.obj calc_dose.obj util.obj
if errorlevel 1 ( echo Build failed at link. & exit /b 1 )

echo.
echo Build succeeded: CCC_DoseEngine_win.exe
exit /b 0
