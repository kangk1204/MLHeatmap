"""Biomarker discovery API endpoints."""

import asyncio
import json
import numpy as np
from fastapi import APIRouter, Query, Request
from fastapi.responses import JSONResponse, StreamingResponse

router = APIRouter(tags=["biomarker"])


@router.get("/biomarker/stream")
async def biomarker_stream(
    request: Request,
    session_id: str = Query(...),
    n_top_genes: int = Query(20, ge=5, le=200),
    n_estimators: int = Query(500, ge=50, le=2000),
    cv_folds: int = Query(5, ge=2, le=10),
    model: str = Query("rf"),
    panel_method: str = Query("forward"),
):
    """Run biomarker analysis with SSE progress streaming."""
    session = request.app.state.sessions.get(session_id)
    if not session or session.normalized is None:
        return JSONResponse({"error": "No normalized data"}, status_code=404)

    def validate_request():
        if len(session.groups) < 2:
            return None, "Need at least 2 groups"

        sample_groups = {}
        for group, samples in session.groups.items():
            for s in samples:
                if s in session.sample_names and s not in session.excluded_samples:
                    sample_groups.setdefault(group, []).append(
                        session.sample_names.index(s)
                    )

        if len(sample_groups) < 2:
            return None, f"Only {len(sample_groups)} group(s) have valid samples after exclusion — need ≥ 2"

        for gname, gidx in sample_groups.items():
            if len(gidx) < 2:
                return None, f"Group '{gname}' has only {len(gidx)} sample(s) after exclusion — need ≥ 2"

        valid_models = {"rf", "random_forest", "xgboost", "xgb", "lightgbm", "lgbm",
                        "logistic_regression", "logistic", "lr_l1",
                        "svm", "svm_linear", "linear_svm"}
        if model.lower().replace(" ", "_") not in valid_models:
            return None, (
                f"Unknown model: {model}. Choose from: rf, xgboost, lightgbm, "
                "logistic_regression, svm_linear"
            )

        valid_panels = {"forward", "lasso", "stability", "mrmr"}
        if panel_method not in valid_panels:
            return None, (
                f"Unknown panel method: {panel_method}. "
                f"Choose from: {', '.join(sorted(valid_panels))}"
            )

        return sample_groups, None

    async def event_generator():
        sample_groups, validation_error = validate_request()
        if validation_error:
            yield f"event: app_error\ndata: {json.dumps({'detail': validation_error})}\n\n"
            return

        progress_queue = asyncio.Queue()
        loop = asyncio.get_running_loop()

        def progress_callback(step, pct, msg):
            loop.call_soon_threadsafe(
                progress_queue.put_nowait,
                {"step": step, "pct": pct, "msg": msg},
            )

        from mlheatmap.core.biomarker import run_biomarker_analysis

        future = loop.run_in_executor(
            None,
            lambda: run_biomarker_analysis(
                expression=session.normalized,
                gene_names=session.gene_names,
                sample_groups=sample_groups,
                n_top_genes=n_top_genes,
                n_estimators=n_estimators,
                cv_folds=cv_folds,
                model=model,
                panel_method=panel_method,
                progress_callback=progress_callback,
            ),
        )

        while not future.done():
            try:
                progress = await asyncio.wait_for(progress_queue.get(), timeout=0.3)
                yield f"event: progress\ndata: {json.dumps(progress)}\n\n"
            except asyncio.TimeoutError:
                continue

        # Drain remaining progress messages after future completes
        while not progress_queue.empty():
            try:
                progress = progress_queue.get_nowait()
                yield f"event: progress\ndata: {json.dumps(progress)}\n\n"
            except asyncio.QueueEmpty:
                break

        try:
            result = future.result()
            session.biomarker_results = result
            yield f"event: complete\ndata: {json.dumps(result, default=_json_safe)}\n\n"
        except Exception as e:
            import traceback
            traceback.print_exc()
            yield f"event: app_error\ndata: {json.dumps({'detail': str(e)})}\n\n"

    return StreamingResponse(event_generator(), media_type="text/event-stream")


@router.get("/biomarker/deg")
async def deg_analysis(
    request: Request,
    session_id: str = Query(...),
    method: str = Query("wilcoxon"),
    log2fc_threshold: float = Query(1.0, ge=0),
    pvalue_threshold: float = Query(0.05, ge=0, le=1),
    use_raw_pvalue: bool = Query(False),
):
    """Run DEG analysis between groups."""
    session = request.app.state.sessions.get(session_id)
    if not session or session.normalized is None:
        return JSONResponse({"error": "No normalized data"}, status_code=404)

    if len(session.groups) != 2:
        return JSONResponse(
            {"error": "DEG requires exactly 2 groups"}, status_code=400
        )

    # Validate method
    valid_methods = ("wilcoxon", "ttest")
    if method not in valid_methods:
        return JSONResponse(
            {"error": f"Unknown DEG method '{method}'. Use one of: {', '.join(valid_methods)}"},
            status_code=400,
        )

    from mlheatmap.core.deg import compute_deg

    # Build sample indices from groups
    sample_groups = {}
    for group, samples in session.groups.items():
        for s in samples:
            if s in session.sample_names and s not in session.excluded_samples:
                sample_groups.setdefault(group, []).append(
                    session.sample_names.index(s)
                )

    # Validate both groups survived exclusion and have enough samples
    if len(sample_groups) != 2:
        missing = [g for g in session.groups if g not in sample_groups]
        return JSONResponse(
            {"error": f"Group(s) {missing} have no valid samples after exclusion"},
            status_code=400,
        )
    for group_name, indices in sample_groups.items():
        if len(indices) < 2:
            return JSONResponse(
                {"error": f"Group '{group_name}' has only {len(indices)} sample(s) after exclusion — need ≥ 2"},
                status_code=400,
            )

    # Use size-factor-normalized counts for log2FC when using DESeq2-VST
    # (raw counts ignore library-size differences; VST is non-linear)
    sf_counts = None
    if session.norm_method == "deseq2" and session.size_factors is not None:
        df = session.mapped_counts if session.mapped_counts is not None else session.raw_counts
        if df is not None:
            raw = df.values.astype(float)
            sf_counts = raw / session.size_factors[np.newaxis, :]

    result = compute_deg(
        expression=session.normalized,
        gene_names=session.gene_names,
        sample_groups=sample_groups,
        method=method,
        log2fc_threshold=log2fc_threshold,
        pvalue_threshold=pvalue_threshold,
        use_raw_pvalue=use_raw_pvalue,
        raw_counts=sf_counts,
    )

    # Store for heatmap use
    session.deg_results = result

    # For ≤20,000 genes, send all results for complete volcano plot
    # For larger datasets, send top by p-value + all significant
    all_results = result["results"]
    if len(all_results) > 20000:
        top_results = all_results[:5000]
        sig_genes = [r for r in all_results if r["direction"] != "ns"]
        seen = {r["gene"] for r in top_results}
        for r in sig_genes:
            if r["gene"] not in seen:
                top_results.append(r)
        all_results = top_results

    response = {
        "results": all_results,
        "summary": result["summary"],
        "group_names": result["group_names"],
        "thresholds": result["thresholds"],
        "method": result["method"],
        "pvalue_type": result.get("pvalue_type", "fdr"),
    }

    return JSONResponse(response, headers={"Content-Type": "application/json"})


def _json_safe(obj):
    """JSON serializer for numpy types."""
    import numpy as np
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")
