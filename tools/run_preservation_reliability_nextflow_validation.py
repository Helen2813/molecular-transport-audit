"""Run preservation reliability through Nextflow and Docker."""

from __future__ import annotations

import hashlib
import json
import shlex
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from run_preservation_nextflow_validation import (
    IMAGE_TAG,
    NEXTFLOW_CONFIG,
    find_wsl,
    to_wsl_path,
    validate_environment,
)


RUNNER_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]

WORKFLOW_FILE = (
    ROOT
    / "workflow"
    / "preservation_reliability_regression.nf"
)

REPORT_DIR = (
    ROOT
    / "reports"
    / "nextflow_preservation_reliability_validation"
)

VERIFICATION_FILE = (
    REPORT_DIR
    / "verification"
    / "reliability_nextflow_verification.json"
)

PROCESS_MANIFEST = (
    REPORT_DIR
    / "preserve_reliability"
    / "reliability_manifest.json"
)

OBSERVED_RELIABILITY = (
    REPORT_DIR
    / "preserve_reliability"
    / "observed_reliability.csv"
)

OBSERVED_LOO = (
    REPORT_DIR
    / "preserve_reliability"
    / "observed_gene_leave_one_out.csv"
)

RUN_MANIFEST = (
    REPORT_DIR
    / "nextflow_preservation_reliability_validation_manifest.json"
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
        "run workflow/"
        "preservation_reliability_regression.nf "
        "-profile docker "
        "-output-dir "
        "reports/"
        "nextflow_preservation_reliability_validation"
    )

    print("")
    print("=" * 80)
    print(
        "Run preservation reliability "
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
            "Nextflow reliability workflow "
            "failed with exit code "
            f"{completed.returncode}."
        )


def load_verification() -> dict:
    if not VERIFICATION_FILE.exists():
        raise FileNotFoundError(
            "Nextflow reliability verification "
            "was not created:\n"
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
            "Nextflow reliability verification "
            "reported passed=false."
        )

    return payload


def main() -> None:
    for path in [
        WORKFLOW_FILE,
        NEXTFLOW_CONFIG,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Required file missing:\n"
                f"{path}"
            )

    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- Nextflow preservation reliability validation"
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

    for path in [
        PROCESS_MANIFEST,
        OBSERVED_RELIABILITY,
        OBSERVED_LOO,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                "Expected published artifact "
                f"was not created:\n{path}"
            )

    checks = verification[
        "checks"
    ]

    run_manifest = {
        "runner_version":
            RUNNER_VERSION,
        "generated_at_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),
        "validation_status":
            "PASS",
        "scientific_component":
            "preservation_reliability",
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
            "summary_metrics_match":
                bool(
                    checks[
                        "summary_metrics_match"
                    ]
                ),
            "valid_repeat_counts_match":
                bool(
                    checks[
                        "valid_repeat_counts_match"
                    ]
                ),
            "common_gene_counts_match":
                bool(
                    checks[
                        "common_gene_counts_match"
                    ]
                ),
            "seed_schedule_match":
                bool(
                    checks[
                        "seed_schedule_match"
                    ]
                ),
            "loo_keys_match":
                bool(
                    checks[
                        "loo_keys_match"
                    ]
                ),
            "loo_correlations_match":
                bool(
                    checks[
                        "loo_correlations_match"
                    ]
                ),
            "maximum_summary_metric_difference":
                float(
                    verification[
                        "maximum_summary_metric_difference"
                    ]
                ),
            "maximum_loo_correlation_difference":
                float(
                    verification[
                        "maximum_loo_correlation_difference"
                    ]
                ),
            "observed_loo_rows":
                int(
                    verification[
                        "observed_loo_rows"
                    ]
                ),
        },
        "published_artifacts": {
            "process_manifest": {
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
            "reliability": {
                "path":
                    str(
                        OBSERVED_RELIABILITY
                        .relative_to(ROOT)
                    ),
                "sha256":
                    sha256_file(
                        OBSERVED_RELIABILITY
                    ),
            },
            "gene_leave_one_out": {
                "path":
                    str(
                        OBSERVED_LOO
                        .relative_to(ROOT)
                    ),
                "sha256":
                    sha256_file(
                        OBSERVED_LOO
                    ),
            },
        },
    }

    RUN_MANIFEST.write_text(
        json.dumps(
            run_manifest,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("")
    print("=" * 80)
    print(
        "Nextflow preservation reliability "
        "validation: PASS"
    )
    print("=" * 80)

    print(
        "Summary metrics: PASS"
    )

    print(
        "Valid repeat counts: PASS"
    )

    print(
        "Common gene counts: PASS"
    )

    print(
        "Seed schedule: PASS"
    )

    print(
        "LOO keys: PASS"
    )

    print(
        "LOO correlations: PASS"
    )

    print(
        "LOO rows: "
        f"{verification['observed_loo_rows']}"
    )

    print(
        "Maximum summary metric "
        "difference: "
        f"{verification['maximum_summary_metric_difference']:.3e}"
    )

    print(
        "Maximum LOO correlation "
        "difference: "
        f"{verification['maximum_loo_correlation_difference']:.3e}"
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
