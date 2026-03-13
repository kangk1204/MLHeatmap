#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
venv_root="$repo_root/.venv"
venv_python="$venv_root/bin/python"
bootstrap_python=0
no_launch=0
selected_python=""
selected_version=""

fail() {
  printf '%s\n' "$1" >&2
  exit 1
}

python_minor_version() {
  local python_cmd="$1"
  "$python_cmd" -c "import sys; print(f'{sys.version_info[0]}.{sys.version_info[1]}')" 2>/dev/null
}

select_python() {
  local candidate version resolved

  for candidate in \
    /opt/homebrew/bin/python3.12 \
    /opt/homebrew/opt/python@3.12/bin/python3.12 \
    /Library/Frameworks/Python.framework/Versions/3.12/bin/python3.12 \
    python3.12 \
    /opt/homebrew/bin/python3.11 \
    /opt/homebrew/opt/python@3.11/bin/python3.11 \
    /Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11 \
    python3.11 \
    python3; do
    if command -v "$candidate" >/dev/null 2>&1; then
      resolved="$(command -v "$candidate")"
    elif [[ -x "$candidate" ]]; then
      resolved="$candidate"
    else
      continue
    fi

    version="$(python_minor_version "$resolved" || true)"
    if [[ "$version" == "3.11" || "$version" == "3.12" ]]; then
      selected_python="$resolved"
      selected_version="$version"
      return 0
    fi
  done

  return 1
}

remove_venv() {
  if [[ -d "$venv_root" ]]; then
    printf 'Removing .venv...\n'
    rm -rf "$venv_root"
  fi
}

ensure_venv() {
  local venv_version

  if [[ -x "$venv_python" ]]; then
    venv_version="$(python_minor_version "$venv_python" || true)"
    if [[ -z "$venv_version" ]]; then
      printf 'Recreating .venv because the existing virtual environment is broken...\n'
      remove_venv
    elif [[ "$venv_version" != "$selected_version" ]]; then
      printf 'Recreating .venv because it uses Python %s and the installer selected Python %s...\n' "$venv_version" "$selected_version"
      remove_venv
    fi
  fi

  if [[ ! -x "$venv_python" ]]; then
    printf 'Creating virtual environment...\n'
    "$selected_python" -m venv "$venv_root" || fail "Failed to create .venv with $selected_python."
  fi
}

upgrade_venv_tooling() {
  printf 'Upgrading pip, setuptools, and wheel...\n'
  "$venv_python" -m pip install --upgrade pip setuptools wheel || fail "Virtual environment tooling upgrade failed."
}

install_application() {
  local attempt=1

  while [[ $attempt -le 2 ]]; do
    ensure_venv
    upgrade_venv_tooling

    printf 'Installing MLHeatmap... (attempt %s/2)\n' "$attempt"
    if "$venv_python" -m pip install -e "$repo_root"; then
      return 0
    fi

    if [[ $attempt -eq 1 ]]; then
      printf 'Installation failed. Recreating .venv and retrying once...\n'
      remove_venv
    fi

    attempt=$((attempt + 1))
  done

  fail "MLHeatmap installation failed after recreating .venv. Close programs that may lock files in .venv and rerun bash ./install-macos.sh."
}

bootstrap_python_install() {
  local brew_cmd=""

  if [[ -x /opt/homebrew/bin/brew ]]; then
    brew_cmd="/opt/homebrew/bin/brew"
  elif command -v brew >/dev/null 2>&1; then
    brew_cmd="$(command -v brew)"
  else
    fail "Bootstrap mode requires Homebrew. Install Homebrew first, then rerun bash ./install-macos.sh --bootstrap-python."
  fi

  printf 'Installing Python 3.12 with Homebrew...\n'
  "$brew_cmd" install python@3.12 || fail "Homebrew failed to install python@3.12."
}

if [[ "$(uname -s)" != "Darwin" ]]; then
  fail "install-macos.sh must be run on macOS."
fi

if [[ "$(uname -m)" != "arm64" ]]; then
  printf 'Warning: Apple Silicon was expected, but detected %s. Continuing anyway.\n' "$(uname -m)"
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bootstrap-python)
      bootstrap_python=1
      ;;
    --no-launch)
      no_launch=1
      ;;
    *)
      fail "Unknown option: $1"
      ;;
  esac
  shift
done

if ! select_python; then
  if [[ $bootstrap_python -eq 0 ]]; then
    fail "Python 3.11 or 3.12 was not found. Install Python 3.12 manually or rerun bash ./install-macos.sh --bootstrap-python. The installer only creates a local .venv and does not install packages into global site-packages."
  fi

  bootstrap_python_install
  select_python || fail "Python bootstrap finished, but no compatible interpreter was detected. Open a new shell and rerun bash ./install-macos.sh."
fi

printf 'Using %s (Python %s)\n' "$selected_python" "$selected_version"
install_application

if [[ $no_launch -eq 1 ]]; then
  printf 'Installation completed. Run bash ./run-macos.sh when you want to start the app.\n'
  exit 0
fi

printf 'Starting MLHeatmap at http://127.0.0.1:8765 ...\n'
exec "$venv_python" -m mlheatmap
