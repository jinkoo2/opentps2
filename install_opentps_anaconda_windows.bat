@echo off

if not defined CONDA_PREFIX (
    echo No anaconda prompt detected.
    echo Run this script in an anaconda prompt.
    exit /b
)

conda env list | findstr /c:"OpenTPS" >nul 2>&1
if %errorlevel% == 1 (
    echo Creating the OpenTPS environment
    call conda create --name OpenTPS python=3.12
)

call conda activate OpenTPS

pip3 install --upgrade pip
pip3 install -e .

echo
echo Start the OpenTPS GUI by running
echo      start_opentps_anaconda_windows.bat