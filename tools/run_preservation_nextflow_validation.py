"""Run deterministic preservation through Nextflow and Docker."""

from __future__ import annotations

import hashlib
import json
import os
import shlex
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


RUNNER_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]

WORKFLOW_FILE = (
    ROOT
    / "workflow"
    / "preservation_regression.nf"
)

NEXTFLOW_CONFIG = (
    ROOT
    / "workflow"
    / "nextflow.config"
)

REPORT_DIR = (
    ROOT
    / "reports"
    / "nextflow_preservation_validation"
)

VERIFICATION_FILE = (
    REPORT_DIR
    / "verification"
    / "preservation_nextflow_verification.json"
)

PROCESS_MANIFEST = (
    REPORT_DIR
    / "preserve_programs"
    / "preservation_manifest.json"
)

RUN_MANIFEST = (
    REPORT_DIR
    / "nextflow_preservation_validation_manifest.json"
)

IMAGE_TAG = (
    "molecular-transport-audit-python:0.1.0"
)


def sha256_file(
    path: Path,
) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def find_wsl() -> str:
    detected = shutil.which(
        "wsl.exe"
    )

    if detected:
        return detected

    system_root = os.environ.get(
        "SystemRoot",
        r"C:\Windows",
    )

    candidate = (
        Path(system_root)
        / "System32"
        / "wsl.exe"
    )

    if candidate.is_file():
        return str(candidate)

    raise RuntimeError(
        "wsl.exe was not found."
    )


def wsl_output(
    wsl: str,
    command: str,
) -> str:
    completed = subprocess.run(
        [
            wsl,
            "-e",
            "bash",
            "-lc",
            command,
        ],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )

    if completed.returncode != 0:
        if completed.stdout:
            print(completed.stdout)

        if completed.stderr:
            print(completed.stderr)

        raise RuntimeError(
            "WSL command failed:\n"
            f"{command}"
        )

    return completed.stdout.strip()


def to_wsl_path(
    wsl: str,
    path: Path,
) -> str:
    windows_path = str(
        path.resolve()
    )

    converted = wsl_output(
        wsl,
        (
            "wslpath -a -u "
            + shlex.quote(
                windows_path
            )
        ),
    )

    if not converted:
        raise RuntimeError(
            "Failed to convert "
            f"path: {windows_path}"
        )

    return converted


def validate_environment(
    wsl: str,
) -> dict[str, str]:
    nextflow_path = wsl_output(
        wsl,
        "command -v nextflow",
    )

    docker_path = wsl_output(
        wsl,
        "command -v docker",
    )

    nextflow_version = wsl_output(
        wsl,
        (
            "nextflow -version "
            "| head -n 5"
        ),
    )

    docker_version = wsl_output(
        wsl,
        (
            "docker version "
            "--format "
            "'{{.Client.Version}}"
            "|{{.Server.Version}}'"
        ),
    )

    image_id = wsl_output(
        wsl,
        (
            "docker image inspect "
            "--format '{{.Id}}' "
            + shlex.quote(
                IMAGE_TAG
            )
        ),
    )

    return {
        "nextflow_path":
            nextflow_path,
        "docker_path":
            docker_path,
        "nextflow_version":
            nextflow_version,
        "docker_version":
            docker_version,
        "image_id":
            image_id,
    }


def run_nextflow(
    wsl: str,
) -> None:
    root_wsl = to_wsl_path(
        wsl,
        ROOT,
    )

    command = (
        "set -euo pipefail; "
        "export NXF_ANSI_LOG=0; "
        f"cd {shlex.quote(root_wsl)}; "
        "nextflow "
        "-c workflow/nextflow.config "
        "run workflow/preservation_regression.nf "
        "-profile docker "
        "-output-dir "
        "reports/nextflow_preservation_validation"
    )

    print("")
    print("=" * 80)
    print(
        "Run preservation regression "
        "through Nextflow"
    )
    print("=" * 80)
    print(
        f"$ {command}"
    )
    print("")

    completed = subprocess.run(
        [
            wsl,
            "-e",
            "bash",
            "-lc",
            command,
        ],
        cwd=ROOT,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "Nextflow preservation "
            "workflow failed with "
            f"exit code "
            f"{completed.returncode}."
        )


def load_verification() -> dict:
    if not VERIFICATION_FILE.exists():
        raise FileNotFoundError(
            "Nextflow preservation "
            "verification was not "
            "created:\n"
            f"{VERIFICATION_FILE}"
        )

    payload = json.loads(
        VERIFICATION_FILE.read_text(
            encoding="utf-8"
        )
    )

    if not bool(
        payload.get(
            "passed",
            False,
        )
    ):
        raise RuntimeError(
            "Nextflow preservation "
            "verification reported "
            "passed=false."
        )

    return payload


def main() -> None:
    for path in [
        WORKFLOW_FILE,
        NEXTFLOW_CONFIG,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Missing required file: "
                f"{path}"
            )

    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- Nextflow preservation validation"
    )
    print("=" * 80)

    print(
        f"Runner version: "
        f"{RUNNER_VERSION}"
    )

    print(
        f"Image: {IMAGE_TAG}"
    )

    wsl = find_wsl()

    environment = (
        validate_environment(
            wsl
        )
    )

    print("")
    print(
        "Validated environment:"
    )

    print(
        "  Nextflow: "
        f"{environment['nextflow_path']}"
    )

    print(
        "  Docker: "
        f"{environment['docker_path']}"
    )

    print(
        "  Docker versions: "
        f"{environment['docker_version']}"
    )

    print(
        "  Image ID: "
        f"{environment['image_id']}"
    )

    if REPORT_DIR.exists():
        shutil.rmtree(
            REPORT_DIR
        )

    REPORT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    run_nextflow(
        wsl
    )

    verification = (
        load_verification()
    )

    if not PROCESS_MANIFEST.exists():
        raise FileNotFoundError(
            "Preservation process manifest "
            "was not published:\n"
            f"{PROCESS_MANIFEST}"
        )

    payload = {
        "runner_version":
            RUNNER_VERSION,
        "generated_at_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),
        "validation_status":
            "PASS",
        "scientific_component":
            (
                "deterministic_direct_"
                "preservation"
            ),
        "workflow": {
            "path":
                str(
                    WORKFLOW_FILE
                    .relative_to(ROOT)
                ),
            "sha256":
                sha256_file(
                    WORKFLOW_FILE
                ),
        },
        "nextflow_config": {
            "path":
                str(
                    NEXTFLOW_CONFIG
                    .relative_to(ROOT)
                ),
            "sha256":
                sha256_file(
                    NEXTFLOW_CONFIG
                ),
        },
        "environment":
            environment,
        "verification": {
            "path":
                str(
                    VERIFICATION_FILE
                    .relative_to(ROOT)
                ),
            "sha256":
                sha256_file(
                    VERIFICATION_FILE
                ),
            "structure_match":
                bool(
                    verification[
                        "structure_match"
                    ]
                ),
            "coverage_match":
                bool(
                    verification[
                        "coverage_match"
                    ]
                ),
            "modules":
                verification[
                    "modules"
                ],
            "max_absolute_metric_difference":
                float(
                    verification[
                        "max_absolute_metric_difference"
                    ]
                ),
        },
        "preservation_process_manifest": {
            "path":
                str(
                    PROCESS_MANIFEST
                    .relative_to(ROOT)
                ),
            "sha256":
                sha256_file(
                    PROCESS_MANIFEST
                ),
        },
    }

    RUN_MANIFEST.write_text(
        json.dumps(
            payload,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("")
    print("=" * 80)
    print(
        "Nextflow preservation "
        "validation: PASS"
    )
    print("=" * 80)

    print(
        "Direct preservation "
        "regression: PASS"
    )

    print(
        "Coverage regression: PASS"
    )

    print(
        "Modules: "
        + ", ".join(
            verification[
                "modules"
            ]
        )
    )

    print(
        "Maximum absolute metric "
        "difference: "
        f"{verification['max_absolute_metric_difference']:.3e}"
    )

    print(
        "Manifest: "
        + str(
            RUN_MANIFEST
            .relative_to(ROOT)
        )
    )

    print("=" * 80)


if __name__ == "__main__":
    main()
