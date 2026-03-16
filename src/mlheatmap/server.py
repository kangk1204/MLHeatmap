"""FastAPI application factory."""

import os
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from mlheatmap import __version__

STATIC_DIR = Path(__file__).parent / "static"


def _cors_origins() -> list[str]:
    configured = os.environ.get("MLHEATMAP_CORS_ORIGINS", "")
    extras = [origin.strip() for origin in configured.split(",") if origin.strip()]
    return extras


LOCAL_ORIGIN_REGEX = r"^https?://(localhost|127\.0\.0\.1|\[::1\])(:\d+)?$"


class _UploadTooLargeError(Exception):
    """Raised when the upload middleware observes a request above the app limit."""


def create_app() -> FastAPI:
    app = FastAPI(title="MLHeatmap", version=__version__)

    app.add_middleware(
        CORSMiddleware,
        allow_origins=_cors_origins(),
        allow_origin_regex=LOCAL_ORIGIN_REGEX,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    @app.middleware("http")
    async def enforce_upload_request_limit(request, call_next):
        from mlheatmap.api.upload import MAX_UPLOAD_BYTES

        if request.url.path != "/api/v1/upload":
            return await call_next(request)

        content_length = request.headers.get("content-length")
        if content_length:
            try:
                if int(content_length) > MAX_UPLOAD_BYTES:
                    return JSONResponse(
                        {"error": f"File too large ({int(content_length) // (1024 * 1024)} MB). Maximum is 500 MB."},
                        status_code=413,
                    )
            except ValueError:
                pass

        original_receive = request._receive
        received = 0

        async def limited_receive():
            nonlocal received
            message = await original_receive()
            if message["type"] == "http.request":
                received += len(message.get("body", b""))
                if received > MAX_UPLOAD_BYTES:
                    raise _UploadTooLargeError(
                        f"File too large ({received // (1024 * 1024)} MB). Maximum is 500 MB."
                    )
            return message

        request._receive = limited_receive
        try:
            return await call_next(request)
        except _UploadTooLargeError as exc:
            return JSONResponse({"error": str(exc)}, status_code=413)

    # Session store
    from mlheatmap.api.session import SessionStore
    app.state.sessions = SessionStore()

    # API routers
    from mlheatmap.api import (
        biomarker,
        capabilities,
        export,
        gene_mapping,
        groups,
        heatmap,
        normalize,
        session_routes,
        upload,
    )

    for router_module in [capabilities, upload, gene_mapping, normalize, heatmap, groups, biomarker, export, session_routes]:
        app.include_router(router_module.router, prefix="/api/v1")

    # Static files
    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

    @app.get("/")
    async def index():
        return FileResponse(STATIC_DIR / "index.html")

    return app
