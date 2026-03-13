"""FastAPI application factory."""

from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from mlheatmap import __version__

STATIC_DIR = Path(__file__).parent / "static"


def create_app() -> FastAPI:
    app = FastAPI(title="MLHeatmap", version=__version__)

    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # Session store
    from mlheatmap.api.session import SessionStore
    app.state.sessions = SessionStore()

    # API routers
    from mlheatmap.api import biomarker, capabilities, export, gene_mapping, groups, heatmap, normalize, upload

    for router_module in [capabilities, upload, gene_mapping, normalize, heatmap, groups, biomarker, export]:
        app.include_router(router_module.router, prefix="/api/v1")

    # Static files
    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

    @app.get("/")
    async def index():
        return FileResponse(STATIC_DIR / "index.html")

    return app
