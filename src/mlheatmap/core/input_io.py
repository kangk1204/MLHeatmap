"""Count-matrix parsing and strict numeric validation helpers."""

from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pandas as pd


class MatrixValidationError(ValueError):
    """Raised when an uploaded matrix is malformed for analysis."""

    def __init__(
        self,
        error: str,
        *,
        invalid_cell_count: int = 0,
        invalid_examples: list[dict[str, str]] | None = None,
        invalid_columns: list[str] | None = None,
    ) -> None:
        super().__init__(error)
        self.error = error
        self.invalid_cell_count = invalid_cell_count
        self.invalid_examples = invalid_examples or []
        self.invalid_columns = invalid_columns or []

    def to_payload(self) -> dict[str, object]:
        payload: dict[str, object] = {"error": self.error}
        if self.invalid_cell_count:
            payload["invalid_cell_count"] = self.invalid_cell_count
            payload["invalid_examples"] = self.invalid_examples
            payload["invalid_columns"] = self.invalid_columns
        return payload


def load_count_matrix_bytes(content: bytes, filename: str) -> pd.DataFrame:
    """Parse a CSV/TSV/XLSX count matrix using the first column as gene ids."""
    name = (filename or "data.csv").lower()
    if name.endswith((".xlsx", ".xls")):
        return pd.read_excel(io.BytesIO(content), index_col=0, engine="openpyxl")
    if name.endswith((".tsv", ".txt")):
        return pd.read_csv(io.BytesIO(content), sep="\t", index_col=0, comment="#")
    return pd.read_csv(io.BytesIO(content), index_col=0, comment="#")


def load_count_matrix_path(path: str | Path) -> pd.DataFrame:
    """Parse a count matrix from disk."""
    file_path = Path(path)
    return load_count_matrix_bytes(file_path.read_bytes(), file_path.name)


def strict_numeric_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Remove fully empty axes and reject any remaining non-numeric content."""
    if df is None:
        raise MatrixValidationError("No data found in file")

    trimmed = df.copy()
    if trimmed.empty:
        raise MatrixValidationError("No valid numeric data found")

    trimmed = trimmed.loc[~trimmed.isna().all(axis=1)]
    trimmed = trimmed.loc[:, ~trimmed.isna().all(axis=0)]
    if trimmed.empty:
        raise MatrixValidationError("No valid numeric data found")

    invalid_examples: list[dict[str, str]] = []
    invalid_columns: set[str] = set()
    invalid_cell_count = 0

    index_series = trimmed.index.to_series()
    missing_gene_mask = index_series.isna() | (index_series.astype(str).str.strip() == "")
    if missing_gene_mask.any():
        missing_rows = index_series[missing_gene_mask]
        invalid_cell_count += int(missing_rows.shape[0])
        invalid_columns.add("gene_id")
        for row_key in missing_rows.index.tolist()[:20]:
            invalid_examples.append(
                {
                    "gene_id": "<missing>",
                    "column": "gene_id",
                    "value": "<missing>",
                }
            )

    numeric = trimmed.apply(pd.to_numeric, errors="coerce")
    invalid_mask = pd.DataFrame(
        ~np.isfinite(numeric.to_numpy(dtype=np.float64)),
        index=trimmed.index,
        columns=trimmed.columns,
    )
    if invalid_mask.to_numpy().any():
        invalid_positions = np.argwhere(invalid_mask.to_numpy())
        invalid_cell_count += int(invalid_positions.shape[0])
        invalid_columns.update(str(trimmed.columns[col_idx]) for _, col_idx in invalid_positions.tolist())
        remaining_slots = max(0, 20 - len(invalid_examples))
        for row_idx, col_idx in invalid_positions[:remaining_slots]:
            raw_value = trimmed.iat[row_idx, col_idx]
            invalid_examples.append(
                {
                    "gene_id": str(trimmed.index[row_idx]),
                    "column": str(trimmed.columns[col_idx]),
                    "value": "<missing>" if pd.isna(raw_value) else str(raw_value),
                }
            )

    if invalid_cell_count:
        raise MatrixValidationError(
            "Input matrix contains missing or non-numeric values. "
            "Remove or fix them before upload.",
            invalid_cell_count=invalid_cell_count,
            invalid_examples=invalid_examples[:20],
            invalid_columns=sorted(invalid_columns),
        )

    return numeric.astype(np.float64)


def filter_low_expression(
    df: pd.DataFrame,
    *,
    min_count: int = 10,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Remove low-expression genes with the app's default threshold."""
    n_before_filter = int(df.shape[0])
    min_samples = max(2, int(df.shape[1] * 0.2))
    expressed_mask = (df >= min_count).sum(axis=1) >= min_samples
    filtered = df.loc[expressed_mask]
    return filtered, {
        "before": n_before_filter,
        "after": int(filtered.shape[0]),
        "removed": n_before_filter - int(filtered.shape[0]),
        "min_count": min_count,
        "min_samples": min_samples,
    }
