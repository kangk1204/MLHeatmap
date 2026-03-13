"""CLI entry point for mlheatmap."""

import argparse
from pathlib import Path
import subprocess
import threading
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


def main():
    parser = argparse.ArgumentParser(
        prog="mlheatmap",
        description="MLHeatmap - Interactive RNA-seq heatmap & biomarker discovery tool",
    )
    parser.add_argument("--port", type=int, default=8765, help="Server port (default: 8765)")
    parser.add_argument("--host", type=str, default="127.0.0.1", help="Server host")
    parser.add_argument("--no-browser", action="store_true", help="Don't auto-open browser")
    args = parser.parse_args()

    url = f"http://{args.host}:{args.port}"

    if not args.no_browser:
        threading.Timer(1.5, _open_browser, args=[url]).start()

    from mlheatmap import __version__
    print(f"\n  MLHeatmap v{__version__}")
    print(f"  Running at {url}")
    print(f"  Press Ctrl+C to stop\n")

    import uvicorn
    from mlheatmap.server import create_app

    app = create_app()
    uvicorn.run(app, host=args.host, port=args.port, log_level="warning")
