"""Upload API endpoints."""

from fastapi import APIRouter, File, Request, UploadFile
from fastapi.responses import JSONResponse

router = APIRouter(tags=["upload"])


@router.post("/upload")
async def upload_file(request: Request, file: UploadFile = File(...)):
    """Upload a count matrix file (CSV/TSV/XLSX)."""
    import io
    import pandas as pd
    from mlheatmap.core.gene_mapping import detect_id_type

    sessions = request.app.state.sessions
    session = sessions.create()

    content = await file.read()
    max_size = 500 * 1024 * 1024  # 500 MB
    if len(content) > max_size:
        return JSONResponse(
            {"error": f"File too large ({len(content) // (1024*1024)} MB). Maximum is 500 MB."},
            status_code=413,
        )
    filename = file.filename or "data.csv"

    try:
        if filename.endswith(".xlsx") or filename.endswith(".xls"):
            df = pd.read_excel(io.BytesIO(content), index_col=0, engine="openpyxl")
        elif filename.endswith(".tsv") or filename.endswith(".txt"):
            df = pd.read_csv(io.BytesIO(content), sep="\t", index_col=0, comment="#")
        else:
            df = pd.read_csv(io.BytesIO(content), index_col=0, comment="#")
    except Exception as e:
        return JSONResponse({"error": f"Failed to parse file: {e}"}, status_code=400)

    # Coerce to numeric
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all", axis=0).dropna(how="all", axis=1)
    df = df.fillna(0)

    # Remove all-zero rows
    df = df.loc[(df != 0).any(axis=1)]

    if df.empty:
        return JSONResponse({"error": "No valid numeric data found"}, status_code=400)

    n_before_filter = df.shape[0]

    # Low-expression gene filtering:
    # Keep genes expressed (count >= 10) in at least 2 samples (or 20% of samples)
    import numpy as np
    min_count = 10
    min_samples = max(2, int(df.shape[1] * 0.2))
    expressed_mask = (df >= min_count).sum(axis=1) >= min_samples
    df = df.loc[expressed_mask]

    n_after_filter = df.shape[0]
    n_filtered_out = n_before_filter - n_after_filter

    gene_ids = df.index.astype(str).tolist()
    species, id_type = detect_id_type(gene_ids)

    session.raw_counts = df
    session.sample_names = df.columns.tolist()
    session.gene_names = gene_ids
    session.species = species
    session.id_type = id_type

    preview = df.head(10).reset_index()
    preview.columns = ["gene_id"] + list(preview.columns[1:])

    return {
        "session_id": session.id,
        "shape": list(df.shape),
        "sample_names": df.columns.tolist(),
        "gene_id_sample": gene_ids[:5],
        "detected_species": species,
        "detected_id_type": id_type,
        "preview": preview.to_dict(orient="records"),
        "filtering": {
            "before": n_before_filter,
            "after": n_after_filter,
            "removed": n_filtered_out,
            "min_count": min_count,
            "min_samples": min_samples,
        },
    }
