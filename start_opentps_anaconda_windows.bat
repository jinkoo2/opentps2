@echo off

if not defined CONDA_PREFIX (
     echo No anaconda prompt detected.
     echo Run this script in an anaconda prompt.
     exit /b
)

REM Check if already in OpenTPS environment, if not try to activate it
if /i not "%CONDA_DEFAULT_ENV%"=="OpenTPS" (
     echo Activating OpenTPS environment...
     call conda activate OpenTPS 2>nul || echo Warning: Could not activate OpenTPS environment
)

set "PYTHONPATH=%~dp0opentps_core;%~dp0opentps_gui"

python "%~dp0opentps_gui\opentps\gui\main.py"