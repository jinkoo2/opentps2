@echo off
REM Build CCC_DoseEngine using Cygwin g++ and make.
REM Requires Cygwin with gcc-g++ and make packages installed.
REM   Run Cygwin Setup, choose "gcc-g++" and "make" under Devel, then install.

set "SCRIPT_DIR=%~dp0"
set "CCC_DIR=%SCRIPT_DIR%"
set "CYGWIN_BIN=C:\cygwin64\bin"

if not exist "%CYGWIN_BIN%\bash.exe" (
  echo ERROR: Cygwin not found at C:\cygwin64\bin. Edit CYGWIN_BIN in this script if needed.
  exit /b 1
)

REM Convert Windows path to Cygwin path (e.g. C:\foo -> /cygdrive/c/foo)
set "CCC_UNIX=%CCC_DIR:\=/%"
set "CCC_UNIX=%CCC_UNIX:C:=/cygdrive/c%"
set "CCC_UNIX=%CCC_UNIX:D:=/cygdrive/d%"
set "CCC_UNIX=%CCC_UNIX:E:=/cygdrive/e%"

echo Building CCC_DoseEngine with Cygwin...
"%CYGWIN_BIN%\bash.exe" -l -c "cd '%CCC_UNIX%' && (make -f Makefile.win 2>/dev/null || (g++ -O2 -c convolution.cpp calc_deff.cpp parse_func.cpp make_poly.cpp raytrace.cpp terma_dose_masks.cpp terma_kerma.cpp calc_dose.cpp util.cpp && g++ -O2 *.o -o CCC_DoseEngine_win.exe && rm -f *.o)) && (test -f CCC_DoseEngine_win && mv -f CCC_DoseEngine_win CCC_DoseEngine_win.exe || true)"
if errorlevel 1 (
  echo.
  echo Build failed. Install Cygwin compiler packages:
  echo   1. Run C:\cygwin64\setup-x86_64.exe
  echo   2. Choose "Install from Internet" or "Download Without Installing"
  echo   3. Select packages: under "Devel" check "gcc-g++" and "make"
  echo   4. Finish setup, then run this script again.
  exit /b 1
)

REM Cygwin g++ may produce CCC_DoseEngine_win (no .exe); ensure .exe for the Python runner
if exist "%CCC_DIR%CCC_DoseEngine_win" if not exist "%CCC_DIR%CCC_DoseEngine_win.exe" (
  move "%CCC_DIR%CCC_DoseEngine_win" "%CCC_DIR%CCC_DoseEngine_win.exe"
)

echo.
echo Build succeeded: CCC_DoseEngine_win.exe
exit /b 0
