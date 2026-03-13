"""CLI entry point for mlheatmap."""

import argparse
import json
from pathlib import Path
import subprocess
import threading
import time
import urllib.error
import urllib.request
import webbrowser


def _is_wsl() -> bool:
    """Detect Windows Subsystem for Linux."""
    proc_version = Path("/proc/version")
    if not proc_version.exists():
        return False
    try:
        return "microsoft" in proc_version.read_text(encoding="utf-8").lower()
    except OSError:
        return False


def _open_browser(url: str) -> None:
    """Open the app URL in a browser on native OSes and WSL."""
    if _is_wsl():
        subprocess.Popen(
            ["cmd.exe", "/c", "start", "", url],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return
    webbrowser.open(url)


def _browser_host(host: str) -> str:
    """Use a local browser-safe host when the server binds to all interfaces."""
    if host in {"0.0.0.0", "::", ""}:
        return "127.0.0.1"
    return host


def _wait_for_server(url: str, timeout_seconds: float = 45.0) -> bool:
    """Poll a local HTTP endpoint until it responds or times out."""
    deadline = time.time() + timeout_seconds
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(url, timeout=2) as response:
                if 200 <= response.status < 500:
                    return True
        except (OSError, urllib.error.URLError):
            time.sleep(0.5)
    return False


def _open_browser_when_ready(app_url: str, readiness_url: str) -> None:
    """Wait for the local server to answer before opening a browser tab."""
    if _wait_for_server(readiness_url):
        _open_browser(app_url)


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="mlheatmap",
        description="MLHeatmap - Interactive RNA-seq heatmap & biomarker discovery tool",
    )
    parser.add_argument("--port", type=int, default=8765, help="Server port (default: 8765)")
    parser.add_argument("--host", type=str, default="127.0.0.1", help="Server host")
    parser.add_argument("--no-browser", action="store_true", help="Don't auto-open browser")
    parser.add_argument(
        "--self-check",
        action="store_true",
        help="Validate installed dependencies and packaged assets, then exit",
    )
    args = parser.parse_args()

    browser_host = _browser_host(args.host)
    url = f"http://{browser_host}:{args.port}"
    readiness_url = f"{url}/api/v1/capabilities"

    if args.self_check:
        from mlheatmap.core.install_check import run_install_self_check

        report = run_install_self_check()
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0

    if not args.no_browser:
        threading.Thread(
            target=_open_browser_when_ready,
            args=(url, readiness_url),
            daemon=True,
        ).start()

    from mlheatmap import __version__
    print(f"\n  MLHeatmap v{__version__}")
    print(f"  Running at {url}")
    print(f"  Press Ctrl+C to stop\n")

    import uvicorn
    from mlheatmap.server import create_app

    app = create_app()
    uvicorn.run(app, host=args.host, port=args.port, log_level="warning")
    return 0
