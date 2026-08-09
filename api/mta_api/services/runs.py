"""Persistent run registry and validated artifact access."""

from __future__ import annotations

import json
import os
import re
import secrets
import subprocess
import sys
from collections import deque
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

from pydantic import ValidationError

from mta_api.config import (
    LEGACY_RUN_ROOT,
    MANAGED_RUNS_ROOT,
    REPOSITORY_ROOT,
    SCIENTIFIC_RUNNER,
    SUPPORTED_ANALYSIS_ID,
)
from mta_api.models import (
    AuditSummary,
    ManagedRunState,
    RunCreateResponse,
    RunDetail,
)


RUN_ID_PATTERN = re.compile(
    r"^[A-Za-z0-9]"
    r"[A-Za-z0-9._-]{0,127}$"
)


class RunArtifactError(RuntimeError):
    """Raised when a persisted run artifact is invalid."""


class RunConflictError(RuntimeError):
    """Raised when another scientific run is active."""


@dataclass(frozen=True)
class RunRecord:
    detail: RunDetail
    summary_path: Path | None
    log_path: Path | None


def utc_now() -> str:
    return datetime.now(
        timezone.utc
    ).isoformat()


def read_json_object(
    path: Path,
) -> dict[str, Any]:
    try:
        payload = json.loads(
            path.read_text(
                encoding="utf-8",
            )
        )
    except (
        OSError,
        json.JSONDecodeError,
    ) as exc:
        raise RunArtifactError(
            f"Could not read JSON artifact: "
            f"{path}"
        ) from exc

    if not isinstance(
        payload,
        dict,
    ):
        raise RunArtifactError(
            f"Expected JSON object: {path}"
        )

    return payload


def atomic_write_json(
    path: Path,
    payload: dict[str, Any],
) -> None:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = path.with_suffix(
        path.suffix + ".tmp"
    )

    temporary.write_text(
        json.dumps(
            payload,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    os.replace(
        temporary,
        path,
    )


def managed_run_root(
    run_id: str,
) -> Path:
    if not RUN_ID_PATTERN.fullmatch(
        run_id
    ):
        raise ValueError(
            "Invalid run ID."
        )

    return (
        MANAGED_RUNS_ROOT
        / run_id
    )


def managed_state_path(
    run_id: str,
) -> Path:
    return (
        managed_run_root(run_id)
        / "run_state.json"
    )


def load_managed_state(
    run_id: str,
) -> ManagedRunState:
    path = managed_state_path(
        run_id
    )

    payload = read_json_object(
        path
    )

    try:
        return (
            ManagedRunState.model_validate(
                payload
            )
        )
    except ValidationError as exc:
        raise RunArtifactError(
            "Invalid managed run state:\n"
            f"{path}\n{exc}"
        ) from exc


def update_managed_state(
    run_id: str,
    **changes: Any,
) -> ManagedRunState:
    state = load_managed_state(
        run_id
    )

    payload = state.model_dump()
    payload.update(
        changes
    )

    try:
        updated = (
            ManagedRunState.model_validate(
                payload
            )
        )
    except ValidationError as exc:
        raise RunArtifactError(
            "Invalid managed run state "
            "transition."
        ) from exc

    atomic_write_json(
        managed_state_path(
            run_id
        ),
        updated.model_dump(),
    )

    return updated


@contextmanager
def registry_lock() -> Iterator[None]:
    MANAGED_RUNS_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    lock_path = (
        MANAGED_RUNS_ROOT
        / ".registry.lock"
    )

    try:
        descriptor = os.open(
            lock_path,
            os.O_CREAT
            | os.O_EXCL
            | os.O_WRONLY,
        )
    except FileExistsError as exc:
        raise RunConflictError(
            "Run registry is currently busy."
        ) from exc

    os.close(
        descriptor
    )

    try:
        yield
    finally:
        lock_path.unlink(
            missing_ok=True
        )


def active_managed_runs() -> list[
    ManagedRunState
]:
    if not MANAGED_RUNS_ROOT.exists():
        return []

    active: list[
        ManagedRunState
    ] = []

    for child in (
        MANAGED_RUNS_ROOT.iterdir()
    ):
        if not child.is_dir():
            continue

        state_path = (
            child
            / "run_state.json"
        )

        if not state_path.exists():
            continue

        state = load_managed_state(
            child.name
        )

        if state.status in {
            "queued",
            "running",
        }:
            active.append(
                state
            )

    return active


def _new_run_id() -> str:
    timestamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%d-%H%M%S"
    )

    suffix = secrets.token_hex(
        3
    )

    return (
        f"gse239948-"
        f"{timestamp}-"
        f"{suffix}"
    )


def create_managed_run(
    analysis_id: str,
) -> RunCreateResponse:
    if (
        analysis_id
        != SUPPORTED_ANALYSIS_ID
    ):
        raise ValueError(
            "Unsupported analysis."
        )

    if not SCIENTIFIC_RUNNER.exists():
        raise RunArtifactError(
            "Scientific runner was not found:\n"
            f"{SCIENTIFIC_RUNNER}"
        )

    with registry_lock():
        active = active_managed_runs()

        if active:
            current = active[0]

            raise RunConflictError(
                "Another scientific run is "
                "already active: "
                f"{current.run_id}"
            )

        run_id = _new_run_id()

        root = managed_run_root(
            run_id
        )

        root.mkdir(
            parents=True,
            exist_ok=False,
        )

        state = ManagedRunState(
            run_id=run_id,
            analysis_id=analysis_id,
            status="queued",
            stage="queued",
            created_at_utc=utc_now(),
            reference_cohort=(
                "GSE238110_DOG2"
            ),
            external_cohort=(
                "GSE239948"
            ),
            outcome_loaded=False,
            pipeline_version="0.1.0",
        )

        atomic_write_json(
            root
            / "run_state.json",
            state.model_dump(),
        )

        bootstrap_path = (
            root
            / "worker_bootstrap.log"
        )

        environment = (
            os.environ.copy()
        )

        environment[
            "PYTHONUNBUFFERED"
        ] = "1"

        command = [
            sys.executable,
            "-m",
            "mta_api.worker",
            "--run-id",
            run_id,
        ]

        popen_options: dict[
            str,
            Any,
        ] = {}

        if os.name == "nt":
            popen_options[
                "creationflags"
            ] = (
                subprocess
                .CREATE_NEW_PROCESS_GROUP
            )
        else:
            popen_options[
                "start_new_session"
            ] = True

        try:
            with bootstrap_path.open(
                "ab"
            ) as handle:
                subprocess.Popen(
                    command,
                    cwd=REPOSITORY_ROOT,
                    stdin=subprocess.DEVNULL,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    env=environment,
                    **popen_options,
                )

        except Exception as exc:
            update_managed_state(
                run_id,
                status="failed",
                stage="failed",
                finished_at_utc=utc_now(),
                error=(
                    "Could not start API worker: "
                    f"{exc}"
                ),
            )

            raise RunArtifactError(
                "Could not launch scientific "
                "worker."
            ) from exc

    return RunCreateResponse(
        run_id=run_id,
        status="queued",
        stage="queued",
        created_at_utc=(
            state.created_at_utc
        ),
    )


def _load_summary(
    path: Path,
) -> AuditSummary:
    payload = read_json_object(
        path
    )

    try:
        return AuditSummary.model_validate(
            payload
        )
    except ValidationError as exc:
        raise RunArtifactError(
            "Summary does not match "
            "API schema contract:\n"
            f"{path}\n{exc}"
        ) from exc


def _legacy_record() -> (
    RunRecord | None
):
    summary_path = (
        LEGACY_RUN_ROOT
        / "summary"
        / "summary.json"
    )

    verification_path = (
        LEGACY_RUN_ROOT
        / "verification"
        / (
            "unified_audit_"
            "verification.json"
        )
    )

    if (
        not summary_path.exists()
        or not verification_path.exists()
    ):
        return None

    verification = read_json_object(
        verification_path
    )

    if (
        verification.get("passed")
        is not True
    ):
        return None

    summary = _load_summary(
        summary_path
    )

    manifest_path = (
        LEGACY_RUN_ROOT
        / (
            "unified_audit_"
            "validation_manifest.json"
        )
    )

    manifest: dict[
        str,
        Any,
    ] = {}

    if manifest_path.exists():
        manifest = read_json_object(
            manifest_path
        )

    pipeline_version_raw = (
        manifest.get(
            "pipeline_version"
        )
    )

    pipeline_version = (
        str(
            pipeline_version_raw
        )
        if pipeline_version_raw
        is not None
        else None
    )

    detail = RunDetail(
        run_id=(
            "unified_audit_validation"
        ),
        analysis_id=(
            summary.run.analysis_id
        ),
        status="completed",
        stage="completed",
        validated=True,
        validation_status="PASS",
        pipeline_version=(
            pipeline_version
        ),
        schema_version=(
            summary.schema_version
        ),
        created_at_utc=(
            summary.generated_at_utc
        ),
        started_at_utc=None,
        finished_at_utc=(
            summary.generated_at_utc
        ),
        generated_at_utc=(
            summary.generated_at_utc
        ),
        reference_cohort=(
            summary.run.reference_cohort
        ),
        external_cohort=(
            summary.run.external_cohort
        ),
        outcome_loaded=(
            summary.run.outcome_loaded
        ),
        module_count=(
            summary.diagnostics.module_count
        ),
        error=None,
    )

    return RunRecord(
        detail=detail,
        summary_path=summary_path,
        log_path=None,
    )


def _managed_record(
    run_id: str,
) -> RunRecord:
    state = load_managed_state(
        run_id
    )

    root = managed_run_root(
        run_id
    )

    summary_path = (
        root
        / "artifacts"
        / "summary"
        / "summary.json"
    )

    summary: (
        AuditSummary | None
    ) = None

    if summary_path.exists():
        summary = _load_summary(
            summary_path
        )
    else:
        summary_path = None

    detail = RunDetail(
        run_id=state.run_id,
        analysis_id=state.analysis_id,
        status=state.status,
        stage=state.stage,
        validated=state.validated,
        validation_status=(
            state.validation_status
        ),
        pipeline_version=(
            state.pipeline_version
        ),
        schema_version=(
            summary.schema_version
            if summary is not None
            else (
                state
                .summary_schema_version
            )
        ),
        created_at_utc=(
            state.created_at_utc
        ),
        started_at_utc=(
            state.started_at_utc
        ),
        finished_at_utc=(
            state.finished_at_utc
        ),
        generated_at_utc=(
            summary.generated_at_utc
            if summary is not None
            else None
        ),
        reference_cohort=(
            state.reference_cohort
        ),
        external_cohort=(
            state.external_cohort
        ),
        outcome_loaded=(
            state.outcome_loaded
        ),
        module_count=(
            summary.diagnostics.module_count
            if summary is not None
            else state.module_count
        ),
        error=state.error,
    )

    log_path = (
        root
        / "run.log"
    )

    if not log_path.exists():
        log_path = None

    return RunRecord(
        detail=detail,
        summary_path=summary_path,
        log_path=log_path,
    )


def get_run_record(
    run_id: str,
) -> RunRecord | None:
    if (
        run_id
        == "unified_audit_validation"
    ):
        return _legacy_record()

    if not RUN_ID_PATTERN.fullmatch(
        run_id
    ):
        return None

    state_path = (
        MANAGED_RUNS_ROOT
        / run_id
        / "run_state.json"
    )

    if not state_path.exists():
        return None

    return _managed_record(
        run_id
    )


def list_run_details() -> list[
    RunDetail
]:
    records: list[
        RunRecord
    ] = []

    legacy = _legacy_record()

    if legacy is not None:
        records.append(
            legacy
        )

    if MANAGED_RUNS_ROOT.exists():
        for child in (
            MANAGED_RUNS_ROOT.iterdir()
        ):
            if (
                not child.is_dir()
                or not (
                    child
                    / "run_state.json"
                ).exists()
            ):
                continue

            records.append(
                _managed_record(
                    child.name
                )
            )

    details = [
        record.detail
        for record in records
    ]

    return sorted(
        details,
        key=lambda item:
        item.created_at_utc,
        reverse=True,
    )


def load_run_summary(
    record: RunRecord,
) -> AuditSummary:
    if record.summary_path is None:
        raise RunArtifactError(
            "Run summary is not available."
        )

    return _load_summary(
        record.summary_path
    )


def tail_run_log(
    record: RunRecord,
    line_count: int,
) -> str:
    if record.log_path is None:
        raise FileNotFoundError(
            "Run log is not available."
        )

    with record.log_path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        lines = deque(
            handle,
            maxlen=line_count,
        )

    return "".join(
        lines
    )
