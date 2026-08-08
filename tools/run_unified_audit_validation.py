"""Run the complete GSE239948 audit through Nextflow and Docker."""

from __future__ import annotations

import hashlib
import json
import shlex
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from run_preservation_docker_validation import (
    build_image,
    image_id,
    resolve_backend,
)

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
    / "molecular_transport_audit_v0_1.nf"
)

REPORT_DIR = (
    ROOT
    / "reports"
    / "unified_audit_validation"
)

SUMMARY_FILE = (
    REPORT_DIR
    / "summary"
    / "summary.json"
)

MODULE_SUMMARY_FILE = (
    REPORT_DIR
    / "summary"
    / "module_summary.csv"
)

VERIFICATION_FILE = (
    REPORT_DIR
    / "verification"
    / "unified_audit_verification.json"
)

RUN_MANIFEST = (
    REPORT_DIR
    / "unified_audit_validation_manifest.json"
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
        "molecular_transport_audit_v0_1.nf "
        "-profile docker "
        "-output-dir "
        "reports/unified_audit_validation"
    )

    print("")
    print("=" * 80)
    print(
        "Run complete molecular "
        "transport audit"
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
            "Unified Nextflow audit failed "
            "with exit code "
            f"{completed.returncode}."
        )


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
        "- unified v0.1 validation"
    )
    print("=" * 80)

    print(
        f"Runner version: "
        f"{RUNNER_VERSION}"
    )

    print("")
    print(
        "Rebuilding scientific Docker image..."
    )

    backend = resolve_backend()

    build_image(
        backend
    )

    current_image_id = image_id(
        backend
    )

    print("")
    print(
        "Scientific image:"
    )
    print(
        f"  Tag: {IMAGE_TAG}"
    )
    print(
        f"  ID: {current_image_id}"
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

    for path in [
        SUMMARY_FILE,
        MODULE_SUMMARY_FILE,
        VERIFICATION_FILE,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                "Expected final artifact "
                f"was not published:\n{path}"
            )

    verification = json.loads(
        VERIFICATION_FILE.read_text(
            encoding="utf-8"
        )
    )

    if not bool(
        verification.get(
            "passed",
            False,
        )
    ):
        raise RuntimeError(
            "Unified audit verifier "
            "reported passed=false."
        )

    summary = json.loads(
        SUMMARY_FILE.read_text(
            encoding="utf-8"
        )
    )

    module_table = pd.read_csv(
        MODULE_SUMMARY_FILE
    )

    run_manifest = {
        "runner_version":
            RUNNER_VERSION,
        "generated_at_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),
        "validation_status":
            "PASS",
        "pipeline_version":
            "0.1.0",
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
        "docker_image": {
            "tag":
                IMAGE_TAG,
            "id":
                current_image_id,
        },
        "environment":
            environment,
        "summary": {
            "path":
                str(
                    SUMMARY_FILE
                    .relative_to(ROOT)
                ),
            "sha256":
                sha256_file(
                    SUMMARY_FILE
                ),
            "schema_version":
                summary[
                    "schema_version"
                ],
        },
        "module_summary": {
            "path":
                str(
                    MODULE_SUMMARY_FILE
                    .relative_to(ROOT)
                ),
            "sha256":
                sha256_file(
                    MODULE_SUMMARY_FILE
                ),
            "rows":
                int(
                    module_table.shape[0]
                ),
        },
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
            "checks":
                verification[
                    "checks"
                ],
            "maximum_reference_metric_difference":
                verification[
                    (
                        "maximum_reference_"
                        "metric_difference"
                    )
                ],
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
        "Unified molecular transport "
        "audit: PASS"
    )
    print("=" * 80)

    print(
        "Reference cohort: "
        f"{summary['run']['reference_cohort']}"
    )

    print(
        "External cohort: "
        f"{summary['run']['external_cohort']}"
    )

    print(
        "Outcome loaded: "
        f"{summary['run']['outcome_loaded']}"
    )

    print("")
    print(
        module_table[
            [
                "module_label",
                "edge_spearman",
                "edge_q_bh_8",
                "loading_spearman",
                "loading_q_bh_8",
                "split_half_median",
                "random_panel_empirical_p",
                "classification_short",
            ]
        ].to_string(
            index=False
        )
    )

    print("")
    print(
        "Summary contract: "
        f"v{summary['schema_version']}"
    )

    print(
        "Summary: "
        + str(
            SUMMARY_FILE
            .relative_to(ROOT)
        )
    )

    print(
        "Verification: "
        + str(
            VERIFICATION_FILE
            .relative_to(ROOT)
        )
    )

    print(
        "Run manifest: "
        + str(
            RUN_MANIFEST
            .relative_to(ROOT)
        )
    )

    print("=" * 80)


if __name__ == "__main__":
    main()
