"""Audit run routes."""

from __future__ import annotations

from fastapi import (
    APIRouter,
    HTTPException,
    Query,
    status,
)
from fastapi.responses import (
    PlainTextResponse,
)

from mta_api.models import (
    AuditSummary,
    RunCreateRequest,
    RunCreateResponse,
    RunDetail,
    RunListResponse,
)
from mta_api.services.runs import (
    RunArtifactError,
    RunConflictError,
    create_managed_run,
    get_run_record,
    list_run_details,
    load_run_summary,
    tail_run_log,
)


router = APIRouter(
    prefix="/runs",
    tags=["runs"],
)


def artifact_error(
    exc: RunArtifactError,
) -> HTTPException:
    return HTTPException(
        status_code=(
            status
            .HTTP_500_INTERNAL_SERVER_ERROR
        ),
        detail={
            "code":
            "artifact_validation_error",
            "message":
            str(exc),
        },
    )


@router.post(
    "",
    response_model=RunCreateResponse,
    status_code=(
        status.HTTP_202_ACCEPTED
    ),
)
def create_run(
    request: RunCreateRequest,
) -> RunCreateResponse:
    try:
        return create_managed_run(
            request.analysis
        )

    except RunConflictError as exc:
        raise HTTPException(
            status_code=(
                status.HTTP_409_CONFLICT
            ),
            detail={
                "code":
                "run_already_active",
                "message":
                str(exc),
            },
        ) from exc

    except RunArtifactError as exc:
        raise artifact_error(
            exc
        ) from exc


@router.get(
    "",
    response_model=RunListResponse,
)
def list_runs() -> RunListResponse:
    try:
        runs = list_run_details()

    except RunArtifactError as exc:
        raise artifact_error(
            exc
        ) from exc

    return RunListResponse(
        runs=runs,
    )


@router.get(
    "/{run_id}",
    response_model=RunDetail,
)
def get_run(
    run_id: str,
) -> RunDetail:
    try:
        record = get_run_record(
            run_id
        )

    except RunArtifactError as exc:
        raise artifact_error(
            exc
        ) from exc

    if record is None:
        raise HTTPException(
            status_code=(
                status.HTTP_404_NOT_FOUND
            ),
            detail="Run not found.",
        )

    return record.detail


@router.get(
    "/{run_id}/summary",
    response_model=AuditSummary,
)
def get_run_summary(
    run_id: str,
) -> AuditSummary:
    try:
        record = get_run_record(
            run_id
        )

    except RunArtifactError as exc:
        raise artifact_error(
            exc
        ) from exc

    if record is None:
        raise HTTPException(
            status_code=(
                status.HTTP_404_NOT_FOUND
            ),
            detail="Run not found.",
        )

    if (
        record.detail.status
        != "completed"
        or not record.detail.validated
        or record.summary_path is None
    ):
        raise HTTPException(
            status_code=(
                status.HTTP_409_CONFLICT
            ),
            detail={
                "code":
                "summary_not_available",
                "message":
                (
                    "Summary becomes available "
                    "only after successful "
                    "verification."
                ),
            },
        )

    try:
        return load_run_summary(
            record
        )

    except RunArtifactError as exc:
        raise artifact_error(
            exc
        ) from exc


@router.get(
    "/{run_id}/log",
    response_class=PlainTextResponse,
)
def get_run_log(
    run_id: str,
    tail: int = Query(
        default=200,
        ge=1,
        le=5000,
    ),
) -> PlainTextResponse:
    try:
        record = get_run_record(
            run_id
        )

    except RunArtifactError as exc:
        raise artifact_error(
            exc
        ) from exc

    if record is None:
        raise HTTPException(
            status_code=(
                status.HTTP_404_NOT_FOUND
            ),
            detail="Run not found.",
        )

    try:
        content = tail_run_log(
            record,
            line_count=tail,
        )

    except FileNotFoundError as exc:
        raise HTTPException(
            status_code=(
                status.HTTP_404_NOT_FOUND
            ),
            detail="Run log not available.",
        ) from exc

    return PlainTextResponse(
        content=content,
    )
