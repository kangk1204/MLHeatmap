@echo off
setlocal

cd /d "%~dp0"

if not exist ".venv\Scripts\python.exe" (
    echo MLHeatmap is not installed in this folder yet.
    echo Run install-windows.cmd first.
    exit /b 1
)

".venv\Scripts\python.exe" -m mlheatmap %*
