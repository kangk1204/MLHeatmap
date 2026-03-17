"""Upload API endpoints."""

from fastapi import APIRouter, File, Request, UploadFile
from fastapi.responses import JSONResponse

router = APIRouter(tags=["upload"])

MAX_UPLOAD_BYTES = 500 * 1024 * 1024


@router.post("/upload")
async def upload_file(request: Request, file: UploadFile = File(...)):
    """Upload a count matrix file (CSV/TSV/XLSX)."""
    from mlheatmap.api.validation import sanitize_upload_filename
    from mlheatmap.core.gene_mapping import detect_id_type
    from mlheatmap.core.input_io import (
        MatrixValidationError,
        load_count_matrix_bytes,
        strict_numeric_matrix,
    )

    content = await file.read(MAX_UPLOAD_BYTES + 1)
    if len(content) > MAX_UPLOAD_BYTES:
        await file.close()
        return JSONResponse(
            {"error": f"File too large ({len(content) // (1024 * 1024)} MB). Maximum is 500 MB."},
            status_code=413,
        )

    filename = sanitize_upload_filename(file.filename)

    try:
        parsed = load_count_matrix_bytes(content, filename)
        df = strict_numeric_matrix(parsed)
    except MatrixValidationError as exc:
        return JSONResponse(exc.to_payload(), status_code=400)
    except Exception as exc:
        return JSONResponse({"error": f"Failed to parse file: {exc}"}, status_code=400)
    finally:
        await file.close()

    n_before_preprocessing = int(df.shape[0])
    df = df.loc[(df != 0).any(axis=1)]
    if df.empty:
        return JSONResponse({"error": "No valid numeric data found"}, status_code=400)
    filtering = {
        "before": n_before_preprocessing,
        "after": int(df.shape[0]),
        "removed": n_before_preprocessing - int(df.shape[0]),
        "applied": False,
        "note": "Low-expression filtering is disabled by default to preserve reproducibility of uploaded matrices.",
    }

    gene_ids = df.index.astype(str).tolist()
    species, id_type = detect_id_type(gene_ids)

    session = request.app.state.sessions.create()
    session.raw_counts = df
    session.sample_names = df.columns.tolist()
    session.gene_names = gene_ids
    session.species = species
    session.id_type = id_type
    session.metadata["upload"] = {
        "filename": filename,
        "shape": [int(df.shape[0]), int(df.shape[1])],
        "detected_species": species,
        "detected_id_type": id_type,
        "filtering": filtering,
        "max_upload_mb": MAX_UPLOAD_BYTES // (1024 * 1024),
    }

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
        "filtering": filtering,
    }
