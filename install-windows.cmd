@echo off
setlocal EnableDelayedExpansion

cd /d "%~dp0"

set "PS_ARGS="

:parse
if "%~1"=="" goto run

if /i "%~1"=="--bootstrap-python" (
    set "PS_ARGS=!PS_ARGS! -BootstrapPython"
) else if /i "%~1"=="--no-launch" (
    set "PS_ARGS=!PS_ARGS! -NoLaunch"
) else (
    set "PS_ARGS=!PS_ARGS! %1"
)

shift
goto parse

:run
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0install-windows.ps1" !PS_ARGS!
