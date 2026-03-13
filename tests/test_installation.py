"""Installation and CLI self-check tests."""

from __future__ import annotations

import json
import subprocess
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer


def test_install_self_check_report():
    from mlheatmap.core.install_check import run_install_self_check

    report = run_install_self_check()

    assert report["capabilities"]["gene_tables"]["human"] is True
    assert report["capabilities"]["gene_tables"]["mouse"] is True
    assert "fastapi" in report["imports"]
    assert "shap" in report["imports"]
    assert "/" in report["routes"]
    assert "/api/v1/capabilities" in report["routes"]
    assert "/static" in report["routes"]


def test_cli_self_check_subprocess():
    proc = subprocess.run(
        [sys.executable, "-m", "mlheatmap", "--self-check"],
        check=True,
        capture_output=True,
        text=True,
    )
    report = json.loads(proc.stdout)

    assert report["capabilities"]["models"]["rf"]["available"] is True
    assert report["capabilities"]["exports"]["results_excel"] is True


def test_browser_host_maps_wildcards_to_loopback():
    from mlheatmap.cli import _browser_host

    assert _browser_host("0.0.0.0") == "127.0.0.1"
    assert _browser_host("::") == "127.0.0.1"
    assert _browser_host("127.0.0.1") == "127.0.0.1"
    assert _browser_host("localhost") == "localhost"


def test_wait_for_server_detects_ready_http_endpoint():
    from mlheatmap.cli import _wait_for_server

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200)
            self.end_headers()
            self.wfile.write(b"ok")

        def log_message(self, format, *args):
            return

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        url = f"http://127.0.0.1:{server.server_address[1]}/health"
        assert _wait_for_server(url, timeout_seconds=3) is True
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
