"""Execute one managed scientific audit run."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import traceback
from pathlib import Path

from mta_api.config import (
    REPOSITORY_ROOT,
    SCIENTIFIC_OUTPUT_ROOT,
    SCIENTIFIC_RUNNER,
)
from mta_api.services.runs import (
    managed_run_root,
    read_json_object,
    update_managed_state,
    utc_now,
)


STAGE_MARKERS = [
    (
        "PRESERVE_PROGRAMS",
        "direct_preservation",
    ),
    (
        "PRESERVE_INFERENCE",
        "permutation_inference",
    ),
    (
        "PRESERVE_RELIABILITY",
        "reliability",
    ),
    (
        "RANDOM_CONTROLS",
        "random_controls",
    ),
    (
        "ASSEMBLE_AUDIT_SUMMARY",
        "summary_assembly",
    ),
    (
        "VERIFY_UNIFIED_AUDIT",
        "verification",
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--run-id",
        required=True,
    )

    return parser.parse_args()


def detect_stage(
    line: str,
) -> str | None:
    for marker, stage in (
        STAGE_MARKERS
    ):
        if marker in line:
            return stage

    return None


def validated_source_artifacts() -> tuple[
    dict,
    dict,
]:
    summary_path = (
        SCIENTIFIC_OUTPUT_ROOT
        / "summary"
        / "summary.json"
    )

    verification_path = (
        SCIENTIFIC_OUTPUT_ROOT
        / "verification"
        / (
            "unified_audit_"
            "verification.json"
        )
    )

    if not summary_path.exists():
        raise RuntimeError(
            "Scientific runner completed "
            "without summary.json."
        )

    if not verification_path.exists():
        raise RuntimeError(
            "Scientific runner completed "
            "without verification artifact."
        )

    summary = read_json_object(
        summary_path
    )

    verification = read_json_object(
        verification_path
    )

    if (
        verification.get("passed")
        is not True
    ):
        raise RuntimeError(
            "Scientific verification "
            "reported passed=false."
        )

    return (
        summary,
        verification,
    )


def copy_artifacts(
    run_root: Path,
) -> None:
    destination = (
        run_root
        / "artifacts"
    )

    if destination.exists():
        shutil.rmtree(
            destination
        )

    shutil.copytree(
        SCIENTIFIC_OUTPUT_ROOT,
        destination,
    )


def execute(
    run_id: str,
) -> None:
    root = managed_run_root(
        run_id
    )

    log_path = (
        root
        / "run.log"
    )

    update_managed_state(
        run_id,
        status="running",
        stage="starting",
        started_at_utc=utc_now(),
        worker_pid=os.getpid(),
        error=None,
    )

    with log_path.open(
        "a",
        encoding="utf-8",
        buffering=1,
    ) as log:
        log.write(
            "Molecular Transport Audit "
            "managed API worker\n"
        )

        log.write(
            f"Run ID: {run_id}\n"
        )

        command = [
            os.fspath(
                Path(
                    os.sys.executable
                )
            ),
            os.fspath(
                SCIENTIFIC_RUNNER
            ),
        ]

        process = subprocess.Popen(
            command,
            cwd=REPOSITORY_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )

        update_managed_state(
            run_id,
            scientific_pid=process.pid,
        )

        current_stage = "starting"

        if process.stdout is None:
            raise RuntimeError(
                "Scientific runner stdout "
                "was not captured."
            )

        for line in process.stdout:
            log.write(
                line
            )

            detected = detect_stage(
                line
            )

            if (
                detected is not None
                and detected
                != current_stage
            ):
                current_stage = detected

                update_managed_state(
                    run_id,
                    stage=detected,
                )

        return_code = (
            process.wait()
        )

        if return_code != 0:
            raise RuntimeError(
                "Scientific runner exited "
                "with code "
                f"{return_code}. "
                "See run.log."
            )

        summary, _ = (
            validated_source_artifacts()
        )

        copy_artifacts(
            root
        )

        module_count = len(
            summary.get(
                "modules",
                [],
            )
        )

        schema_version = str(
            summary.get(
                "schema_version",
                "",
            )
        )

        if (
            schema_version
            != "0.1.0"
        ):
            raise RuntimeError(
                "Unexpected summary schema: "
                f"{schema_version}"
            )

        update_managed_state(
            run_id,
            status="completed",
            stage="completed",
            finished_at_utc=utc_now(),
            validated=True,
            validation_status="PASS",
            summary_schema_version=(
                schema_version
            ),
            module_count=(
                module_count
            ),
            error=None,
        )

        log.write(
            "\nManaged API run: PASS\n"
        )


def main() -> None:
    args = parse_args()

    try:
        execute(
            args.run_id
        )

    except Exception as exc:
        root = managed_run_root(
            args.run_id
        )

        log_path = (
            root
            / "run.log"
        )

        with log_path.open(
            "a",
            encoding="utf-8",
        ) as log:
            log.write(
                "\nManaged API run: FAIL\n"
            )

            traceback.print_exc(
                file=log
            )

        try:
            update_managed_state(
                args.run_id,
                status="failed",
                stage="failed",
                finished_at_utc=utc_now(),
                validated=False,
                validation_status=None,
                error=str(exc),
            )
        except Exception:
            pass

        raise


if __name__ == "__main__":
    main()
