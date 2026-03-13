"""Runtime capability inspection for optional features and packaged data."""

from __future__ import annotations

import importlib
import platform
from functools import lru_cache
from typing import Any

from mlheatmap.core.gene_mapping import has_gene_table

SUPPORTED_PYTHON = ">=3.11,<3.13"

MODEL_SPECS: dict[str, dict[str, str | None]] = {
    "rf": {"label": "Random Forest", "module": None, "install_profile": "core"},
    "logistic_regression": {
        "label": "Logistic Regression (L1)",
        "module": None,
        "install_profile": "core",
    },
    "svm_linear": {"label": "SVM (Linear)", "module": None, "install_profile": "core"},
    "xgboost": {"label": "XGBoost", "module": "xgboost", "install_profile": "full"},
    "lightgbm": {"label": "LightGBM", "module": "lightgbm", "install_profile": "full"},
}

MODEL_ALIASES = {
    "rf": "rf",
    "random_forest": "rf",
    "logistic_regression": "logistic_regression",
    "logistic": "logistic_regression",
    "lr_l1": "logistic_regression",
    "svm": "svm_linear",
    "svm_linear": "svm_linear",
    "linear_svm": "svm_linear",
    "xgboost": "xgboost",
    "xgb": "xgboost",
    "lightgbm": "lightgbm",
    "lgbm": "lightgbm",
}


def normalize_model_name(model_name: str | None) -> str:
    """Normalize model aliases to canonical public ids."""
    key = (model_name or "rf").lower().replace(" ", "_")
    return MODEL_ALIASES.get(key, key)


@lru_cache(maxsize=None)
def _module_available(module_name: str) -> tuple[bool, str | None]:
    try:
        importlib.import_module(module_name)
        return True, None
    except Exception as exc:
        return False, str(exc)


def get_model_capability(model_name: str) -> dict[str, Any]:
    """Return availability metadata for a model id or alias."""
    model_id = normalize_model_name(model_name)
    spec = MODEL_SPECS.get(model_id)
    if spec is None:
        return {
            "id": model_id,
            "label": model_id,
            "available": False,
            "known": False,
            "install_profile": None,
            "unavailable_reason": f"Unknown model: {model_name}",
        }

    module_name = spec["module"]
    available = True
    if module_name:
        available, _ = _module_available(module_name)

    reason = None
    if not available:
        reason = (
            f"{spec['label']} is not available in this installation. "
            "Install mlheatmap[full] or use Docker for the full model set."
        )

    return {
        "id": model_id,
        "label": spec["label"],
        "available": available,
        "known": True,
        "install_profile": spec["install_profile"],
        "unavailable_reason": reason,
    }


def get_capabilities() -> dict[str, Any]:
    """Return runtime capabilities for the current installation."""
    models = {model_id: get_model_capability(model_id) for model_id in MODEL_SPECS}
    return {
        "python": {
            "version": platform.python_version(),
            "supported": SUPPORTED_PYTHON,
        },
        "gene_tables": {
            "human": has_gene_table("human"),
            "mouse": has_gene_table("mouse"),
        },
        "models": models,
        "available_models": [model_id for model_id, info in models.items() if info["available"]],
        "exports": {
            "image_mode": "browser",
            "results_excel": True,
            "server_image_exports": False,
            "formats": ["png", "svg"],
        },
        "docker": {
            "recommended_for": ["xgboost", "lightgbm"],
            "entrypoint": ["mlheatmap", "--host", "0.0.0.0", "--no-browser"],
        },
    }
