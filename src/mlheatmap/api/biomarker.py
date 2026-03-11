"""Biomarker discovery API endpoints."""

import asyncio
import json
from fastapi import APIRouter, Query, Request
from fastapi.responses import JSONResponse, StreamingResponse

router = APIRouter(tags=["biomarker"])


@router.get("/biomarker/stream")
async def biomarker_stream(
    request: Request,
    session_id: str = Query(...),
    n_top_genes: int = Query(20, ge=5, le=100),
    n_estimators: int = Query(500, ge=50, le=2000),
    cv_folds: int = Query(5, ge=2, le=10),
    model: str = Query("rf"),
):
    """Run biomarker analysis with SSE progress streaming."""
    session = request.app.state.sessions.get(session_id)
    if not session or session.normalized is None:
        return JSONResponse({"error": "No normalized data"}, status_code=404)

    if len(session.groups) < 2:
        return JSONResponse({"error": "Need at least 2 groups"}, status_code=400)

    async def event_generator():
        progress_queue = asyncio.Queue()
        loop = asyncio.get_running_loop()

        def progress_callback(step, pct, msg):
            loop.call_soon_threadsafe(
                progress_queue.put_nowait,
                {"step": step, "pct": pct, "msg": msg},
            )

        from mlheatmap.core.biomarker import run_biomarker_analysis

        # Build sample indices from groups
        sample_groups = {}
        for group, samples in session.groups.items():
            for s in samples:
                if s in session.sample_names:
                    sample_groups.setdefault(group, []).append(
                        session.sample_names.index(s)
                    )

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
            yield f"event: error\ndata: {json.dumps({'detail': str(e)})}\n\n"

    return StreamingResponse(event_generator(), media_type="text/event-stream")


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
