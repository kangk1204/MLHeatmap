param(
    [switch]$NoLaunch,
    [switch]$BootstrapPython
)

$ErrorActionPreference = "Stop"
$repoRoot = $PSScriptRoot
$venvPython = Join-Path $repoRoot ".venv\Scripts\python.exe"
$venvRoot = Join-Path $repoRoot ".venv"

function Fail {
    param(
        [string]$Message
    )

    Write-Host $Message -ForegroundColor Red
    exit 1
}

function Invoke-CommandArray {
    param(
        [string[]]$Command,
        [string[]]$Arguments = @()
    )

    $exe = $Command[0]
    $baseArgs = @()
    if ($Command.Length -gt 1) {
        $baseArgs = $Command[1..($Command.Length - 1)]
    }

    & $exe @baseArgs @Arguments
    return $LASTEXITCODE
}

function Remove-Venv {
    if (Test-Path $venvRoot) {
        Write-Host "Removing .venv..."
        Remove-Item $venvRoot -Recurse -Force
    }
}

function Ensure-Venv {
    param(
        [object]$Python
    )

    if (Test-Path $venvPython) {
        $venvVersion = Get-PythonMinorVersion -Command @($venvPython)
        if (-not $venvVersion) {
            Write-Host "Recreating .venv because the existing virtual environment is broken..."
            Remove-Venv
        } elseif ($venvVersion -ne $Python.Version) {
            Write-Host "Recreating .venv because it uses Python $venvVersion and the installer selected Python $($Python.Version)..."
            Remove-Venv
        }
    }

    if (-not (Test-Path $venvPython)) {
        Write-Host "Creating virtual environment..."
        $exitCode = Invoke-CommandArray -Command $Python.Command -Arguments @("-m", "venv", ".venv")
        if ($exitCode -ne 0) {
            Fail "Failed to create .venv with $($Python.Label)."
        }
    }
}

function Upgrade-VenvTooling {
    Write-Host "Upgrading pip, setuptools, and wheel..."
    & $venvPython -m pip install --upgrade pip setuptools wheel
    if ($LASTEXITCODE -ne 0) {
        Fail "Virtual environment tooling upgrade failed."
    }
}

function Install-Application {
    param(
        [object]$Python
    )

    $attempt = 1
    while ($attempt -le 2) {
        Ensure-Venv -Python $Python
        Upgrade-VenvTooling

        Write-Host "Installing MLHeatmap... (attempt $attempt/2)"
        & $venvPython -m pip install $repoRoot
        if ($LASTEXITCODE -eq 0) {
            return
        }

        if ($attempt -eq 1) {
            Write-Host "Installation failed. Recreating .venv and retrying once..."
            Remove-Venv
        }

        $attempt += 1
    }

    Fail "MLHeatmap installation failed after recreating .venv. Close programs that may lock files in .venv and rerun install-windows.cmd."
}

function Invoke-SelfCheck {
    Write-Host "Running install self-check..."
    & $venvPython -m mlheatmap --self-check
    if ($LASTEXITCODE -ne 0) {
        Fail "MLHeatmap self-check failed."
    }
}

function Get-PythonMinorVersion {
    param(
        [string[]]$Command
    )

    try {
        $exe = $Command[0]
        $baseArgs = @()
        if ($Command.Length -gt 1) {
            $baseArgs = $Command[1..($Command.Length - 1)]
        }

        $output = & $exe @baseArgs -c "import sys; print(str(sys.version_info[0]) + '.' + str(sys.version_info[1]))" 2>$null
        $exitCode = $LASTEXITCODE
        if ($exitCode -eq 0) {
            return ($output | Select-Object -First 1).Trim()
        }
    } catch {
    }

    return $null
}

function Get-CompatiblePython {
    $candidates = @()

    if (Get-Command py -ErrorAction SilentlyContinue) {
        $candidates += [pscustomobject]@{ Command = @("py", "-3.12"); Label = "py -3.12" }
        $candidates += [pscustomobject]@{ Command = @("py", "-3.11"); Label = "py -3.11" }
    }

    $standardPaths = @(
        (Join-Path $env:LocalAppData "Programs\Python\Python312\python.exe"),
        (Join-Path $env:LocalAppData "Programs\Python\Python311\python.exe"),
        (Join-Path $env:ProgramFiles "Python312\python.exe"),
        (Join-Path $env:ProgramFiles "Python311\python.exe")
    )

    foreach ($path in $standardPaths) {
        if ($path -and (Test-Path $path)) {
            $candidates += [pscustomobject]@{ Command = @($path); Label = $path }
        }
    }

    if (Get-Command python -ErrorAction SilentlyContinue) {
        $candidates += [pscustomobject]@{ Command = @("python"); Label = "python" }
    }

    foreach ($candidate in $candidates) {
        $version = Get-PythonMinorVersion -Command $candidate.Command
        if ($version -in @("3.11", "3.12")) {
            return [pscustomobject]@{
                Command = $candidate.Command
                Label = $candidate.Label
                Version = $version
            }
        }
    }

    return $null
}

function Install-BootstrapPython {
    $winget = Get-Command winget -ErrorAction SilentlyContinue
    if (-not $winget) {
        Fail "A compatible Python interpreter was not found and winget is unavailable. Install Python 3.12 manually from https://www.python.org/downloads/windows/ and rerun install-windows.cmd."
    }

    Write-Host "Installing Python 3.12 for the current user with winget..."
    Write-Host "This installs Python itself in the user profile, but packages will still be installed only in .venv."

    & $winget.Source install --id Python.Python.3.12 --exact --scope user --silent --disable-interactivity --accept-package-agreements --accept-source-agreements --custom "InstallAllUsers=0 PrependPath=0 Include_test=0"
    if ($LASTEXITCODE -ne 0) {
        Fail "Python 3.12 bootstrap installation failed."
    }

    Start-Sleep -Seconds 2
}

function Wait-HttpReady {
    param(
        [string]$Url,
        [int]$TimeoutSeconds = 45
    )

    $deadline = (Get-Date).AddSeconds($TimeoutSeconds)
    while ((Get-Date) -lt $deadline) {
        try {
            $response = Invoke-WebRequest -Uri $Url -UseBasicParsing -TimeoutSec 3
            if ($response.StatusCode -ge 200 -and $response.StatusCode -lt 500) {
                return $true
            }
        } catch {
        }

        Start-Sleep -Seconds 1
    }

    return $false
}

$python = Get-CompatiblePython
if (-not $python) {
    Write-Host "Python 3.11 or 3.12 was not found. Attempting bootstrap install..."
    Install-BootstrapPython
    $python = Get-CompatiblePython
    if (-not $python) {
        Fail "Python 3.12 bootstrap installation finished, but no compatible interpreter was detected. Open a new terminal and rerun install-windows.cmd."
    }
}

Write-Host "Using $($python.Label) (Python $($python.Version))"
Install-Application -Python $python
Invoke-SelfCheck

if ($NoLaunch) {
    Write-Host "Installation completed. Run .\run-windows.cmd when you want to start the app."
    exit 0
}

Write-Host "Starting MLHeatmap at http://127.0.0.1:8765 ..."
Start-Process -FilePath $venvPython -ArgumentList @("-m", "mlheatmap", "--no-browser") -WorkingDirectory $repoRoot | Out-Null
if (-not (Wait-HttpReady -Url "http://127.0.0.1:8765/api/v1/capabilities")) {
    Fail "MLHeatmap started, but the local server did not become ready at http://127.0.0.1:8765 within 45 seconds."
}
Start-Process "http://127.0.0.1:8765"
