#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
venv_python="$repo_root/.venv/bin/python"

if [[ ! -x "$venv_python" ]]; then
  printf 'MLHeatmap is not installed in this folder yet.\nRun bash ./install-macos.sh first.\n' >&2
  exit 1
fi

exec "$venv_python" -m mlheatmap "$@"
