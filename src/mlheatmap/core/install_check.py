"""Installation self-checks for packaged assets and required modules."""

from __future__ import annotations

import importlib
import platform
from importlib import resources
from typing import Any

from mlheatmap.core.capabilities import get_capabilities
from mlheatmap.server import create_app

REQUIRED_IMPORTS: dict[str, tuple[str, ...]] = {
    "fastapi": ("fastapi",),
    "uvicorn": ("uvicorn",),
    "numpy": ("numpy",),
    "scipy": ("scipy",),
    "pandas": ("pandas",),
    "scikit-learn": ("sklearn",),
    "shap": ("shap",),
    "openpyxl": ("openpyxl",),
    "python-multipart": ("python_multipart", "multipart"),
    "matplotlib": ("matplotlib",),
    "seaborn": ("seaborn",),
}

STATIC_ASSETS = (
    "index.html",
    "css/main.css",
    "js/api.js",
    "js/app.js",
    "js/biomarker.js",
    "js/export.js",
    "js/groups.js",
    "js/heatmap.js",
    "js/upload.js",
)


def run_install_self_check() -> dict[str, Any]:
    """Validate the installed application without starting the server."""
    import_errors: dict[str, str] = {}
    imported_modules: list[str] = []
    for label, module_names in REQUIRED_IMPORTS.items():
        last_error = None
        for module_name in module_names:
            try:
                importlib.import_module(module_name)
                imported_modules.append(label)
                break
            except Exception as exc:  # pragma: no cover - exercised via failure path in installers
                last_error = exc
        else:
            import_errors[label] = str(last_error)

    if import_errors:
        formatted = ", ".join(f"{label}: {error}" for label, error in sorted(import_errors.items()))
        raise RuntimeError(f"Required Python modules failed to import: {formatted}")

    static_root = resources.files("mlheatmap").joinpath("static")
    missing_assets = [asset for asset in STATIC_ASSETS if not static_root.joinpath(asset).is_file()]
    if missing_assets:
        raise RuntimeError(
            "Packaged static assets are missing: " + ", ".join(missing_assets)
        )

    app = create_app()
    route_paths = {getattr(route, "path", "") for route in app.routes}
    for required_path in ("/", "/api/v1/capabilities", "/static"):
        if required_path not in route_paths:
            raise RuntimeError(f"Application route is missing: {required_path}")

    capabilities = get_capabilities()
    if not capabilities["gene_tables"]["human"] or not capabilities["gene_tables"]["mouse"]:
        raise RuntimeError("Packaged gene mapping tables are missing.")
    if not capabilities["models"]["rf"]["available"]:
        raise RuntimeError("Core Random Forest model is not available.")

    return {
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "imports": imported_modules,
        "static_assets": list(STATIC_ASSETS),
        "routes": sorted(route_paths),
        "capabilities": capabilities,
    }
