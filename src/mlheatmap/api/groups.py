"""Group assignment API endpoints."""

from __future__ import annotations

from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from mlheatmap.api.validation import acquire_session_lease_or_error, get_session_or_error

router = APIRouter(tags=["groups"])


class GroupsRequest(BaseModel):
    session_id: str
    groups: dict


class ExcludeRequest(BaseModel):
    session_id: str
    samples: list


@router.post("/groups")
async def set_groups(request: Request, body: GroupsRequest):
    """Set sample group assignments."""
    lease, error, _ = acquire_session_lease_or_error(request, body.session_id)
    if error is not None:
        return error
    session = lease.session

    try:
        with session.state_lock:
            valid_samples = set(session.sample_names) if session.sample_names else set()

        seen: dict[str, str] = {}
        sanitized_groups: dict[str, list[str]] = {}
        for group_name, samples in body.groups.items():
            unique_in_group = []
            for sample in samples:
                if valid_samples and sample not in valid_samples:
                    return JSONResponse(
                        {"error": f"Unknown sample '{sample}' in group '{group_name}' - not in the uploaded matrix"},
                        status_code=400,
                    )
                if sample in seen:
                    return JSONResponse(
                        {"error": f"Sample '{sample}' appears in both '{seen[sample]}' and '{group_name}'"},
                        status_code=400,
                    )
                seen[sample] = group_name
                if sample not in unique_in_group:
                    unique_in_group.append(sample)
            sanitized_groups[group_name] = unique_in_group

        with session.state_lock:
            session.groups = sanitized_groups
            session.excluded_group_assignments = {}
            session._invalidate_analysis_unlocked()
            groups = dict(session.groups)

        return {"groups": groups, "n_groups": len(groups)}
    finally:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)


@router.get("/groups")
async def get_groups(request: Request, session_id: str):
    """Get current group assignments."""
    session, error, _ = get_session_or_error(request, session_id)
    if error is not None:
        return error
    with session.state_lock:
        return {"groups": dict(session.groups)}


@router.post("/groups/exclude")
async def exclude_samples(request: Request, body: ExcludeRequest):
    """Exclude samples from analysis."""
    lease, error, _ = acquire_session_lease_or_error(request, body.session_id)
    if error is not None:
        return error
    session = lease.session

    try:
        with session.state_lock:
            sample_to_group = {}
            for group_name, group_samples in session.groups.items():
                for sample in group_samples:
                    sample_to_group[sample] = group_name
            for sample in body.samples:
                if sample in sample_to_group:
                    session.excluded_group_assignments[sample] = sample_to_group[sample]
            session.excluded_samples = list(set(session.excluded_samples + body.samples))
            for group_name in session.groups:
                session.groups[group_name] = [
                    sample for sample in session.groups[group_name] if sample not in body.samples
                ]
            session._invalidate_analysis_unlocked()
            excluded = list(session.excluded_samples)
            remaining = [sample for sample in session.sample_names if sample not in session.excluded_samples]

        return {"excluded": excluded, "remaining": remaining}
    finally:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)


@router.post("/groups/include")
async def include_samples(request: Request, body: ExcludeRequest):
    """Re-include previously excluded samples."""
    lease, error, _ = acquire_session_lease_or_error(request, body.session_id)
    if error is not None:
        return error
    session = lease.session

    try:
        with session.state_lock:
            session.excluded_samples = [sample for sample in session.excluded_samples if sample not in body.samples]
            current_assignments = {
                sample: group_name
                for group_name, group_samples in session.groups.items()
                for sample in group_samples
            }
            restored_groups: dict[str, str] = {}
            unassigned_samples: list[str] = []
            for sample in body.samples:
                original_group = session.excluded_group_assignments.pop(sample, None)
                if original_group and original_group in session.groups and sample not in current_assignments:
                    session.groups[original_group].append(sample)
                    restored_groups[sample] = original_group
                elif sample not in current_assignments:
                    unassigned_samples.append(sample)
            session._invalidate_analysis_unlocked()
            excluded = list(session.excluded_samples)
            remaining = [sample for sample in session.sample_names if sample not in session.excluded_samples]

        return {
            "excluded": excluded,
            "remaining": remaining,
            "restored_groups": restored_groups,
            "unassigned_samples": unassigned_samples,
        }
    finally:
        request.app.state.sessions.end_use(lease.session_id, lease.operation_id)
