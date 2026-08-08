"""Validate preservation reliability inside the Docker environment."""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

from run_preservation_docker_validation import (
    IMAGE_TAG,
    DockerBackend,
    build_image,
    container_versions,
    docker_versions,
    image_id,
    bind_mount,
    resolve_backend,
    run_docker,
)


RUNNER_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]

TOOLS_DIR = ROOT / "tools"
TESTS_DIR = ROOT / "tests"

FIXTURE_DIR = (
    ROOT
    / "reference_results"
    / "osteosarcoma_locked"
    / "preservation_fixture"
)

REPORT_ROOT = (
    ROOT
    / "reports"
    / "docker_preservation_reliability_validation"
)

VERIFICATION_FILE = (
    REPORT_ROOT
    / "preservation_reliability_regression"
    / "GSE239948_reliability_verification.json"
)

OUTPUT_MANIFEST = (
    REPORT_ROOT
    / "docker_preservation_reliability_validation_manifest.json"
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


def require_path(
    path: Path,
    description: str,
) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"{description} not found:\n"
            f"{path}"
        )


def run_container_validation(
    backend: DockerBackend,
) -> None:
    print("")
    print("=" * 80)
    print(
        "Run preservation reliability "
        "suite inside Docker"
    )
    print("=" * 80)

    if REPORT_ROOT.exists():
        shutil.rmtree(
            REPORT_ROOT
        )

    REPORT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    command = [
        "run",
        "--rm",
        "--network",
        "none",
        "--read-only",
        "--tmpfs",
        (
            "/tmp:"
            "rw,nosuid,nodev,"
            "size=64m"
        ),
        "--mount",
        bind_mount(
            backend,
            TOOLS_DIR,
            "/opt/mta/tools",
            readonly=True,
        ),
        "--mount",
        bind_mount(
            backend,
            TESTS_DIR,
            "/opt/mta/tests",
            readonly=True,
        ),
        "--mount",
        bind_mount(
            backend,
            FIXTURE_DIR,
            (
                "/opt/mta/"
                "reference_results/"
                "osteosarcoma_locked/"
                "preservation_fixture"
            ),
            readonly=True,
        ),
        "--mount",
        bind_mount(
            backend,
            REPORT_ROOT,
            "/opt/mta/reports",
            readonly=False,
        ),
        IMAGE_TAG,
        "python3",
        (
            "tools/"
            "run_preservation_reliability_validation_suite.py"
        ),
    ]

    run_docker(
        backend,
        command,
    )


def load_verification() -> dict:
    if not VERIFICATION_FILE.exists():
        raise FileNotFoundError(
            "Docker reliability verification "
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
            "Docker preservation reliability "
            "verification reported "
            "passed=false."
        )

    return payload


def write_manifest(
    *,
    backend: DockerBackend,
    docker_info: dict[str, str],
    image_sha: str,
    package_versions: dict[str, str],
    verification: dict,
) -> None:
    checks = verification[
        "checks"
    ]

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
            "preservation_reliability",
        "docker_backend": {
            "kind":
                backend.kind,
            "description":
                backend.description,
            "launcher":
                list(
                    backend.launcher
                ),
        },
        "image": {
            "tag":
                IMAGE_TAG,
            "id":
                image_sha,
        },
        "docker":
            docker_info,
        "container_environment":
            package_versions,
        "isolation": {
            "network":
                "none",
            "root_filesystem":
                "read_only",
            "tools_mount":
                "read_only",
            "tests_mount":
                "read_only",
            "fixture_mount":
                "read_only",
            "report_mount":
                "read_write",
        },
        "reliability_contract": {
            "split_half_repeats":
                int(
                    verification[
                        "split_half_repeats"
                    ]
                ),
            "outcome_loaded":
                bool(
                    verification[
                        "outcome_loaded"
                    ]
                ),
            "expected_seeds": {
                "M34": 44,
                "M11": 144,
                "M24": 244,
                "M40": 344,
            },
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
            "fixture_hashes_match":
                bool(
                    verification[
                        "fixture_hashes_match"
                    ]
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
    }

    OUTPUT_MANIFEST.write_text(
        json.dumps(
            payload,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    require_path(
        (
            TOOLS_DIR
            / "run_preservation_reliability_validation_suite.py"
        ),
        "Reliability validation suite",
    )

    require_path(
        (
            TOOLS_DIR
            / "run_preservation_reliability_regression.py"
        ),
        "Reliability regression runner",
    )

    require_path(
        (
            TESTS_DIR
            / "unit"
            / "test_preservation_reliability.py"
        ),
        "Reliability unit tests",
    )

    require_path(
        (
            FIXTURE_DIR
            / "reliability_fixture_manifest.json"
        ),
        "Reliability fixture manifest",
    )

    require_path(
        (
            FIXTURE_DIR
            / "expected_gene_leave_one_out.csv"
        ),
        "Gene leave-one-out fixture",
    )

    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- Docker preservation reliability validation"
    )
    print("=" * 80)

    print(
        f"Runner version: "
        f"{RUNNER_VERSION}"
    )

    print(
        f"Image: "
        f"{IMAGE_TAG}"
    )

    backend = resolve_backend()

    print("")
    print(
        "Selected Docker backend:"
    )
    print(
        f"  Kind: "
        f"{backend.kind}"
    )
    print(
        f"  {backend.description}"
    )

    docker_info = (
        docker_versions(
            backend
        )
    )

    print(
        "Docker client: "
        f"{docker_info['client']}"
    )
    print(
        "Docker server: "
        f"{docker_info['server']}"
    )

    # Rebuild because preservation_reliability.py
    # is a new scientific-core module.
    build_image(
        backend
    )

    current_image_id = (
        image_id(
            backend
        )
    )

    versions = (
        container_versions(
            backend
        )
    )

    print("")
    print(
        "Container environment:"
    )

    for key, value in (
        versions.items()
    ):
        print(
            f"  {key}: {value}"
        )

    run_container_validation(
        backend
    )

    verification = (
        load_verification()
    )

    write_manifest(
        backend=backend,
        docker_info=docker_info,
        image_sha=current_image_id,
        package_versions=versions,
        verification=verification,
    )

    checks = verification[
        "checks"
    ]

    print("")
    print("=" * 80)
    print(
        "Docker preservation reliability "
        "validation: PASS"
    )
    print("=" * 80)

    print(
        "Fixture hashes: PASS"
    )
    print(
        "Reliability unit tests: PASS"
    )

    print(
        "Summary metrics: "
        + (
            "PASS"
            if checks[
                "summary_metrics_match"
            ]
            else "FAIL"
        )
    )

    print(
        "Valid repeat counts: "
        + (
            "PASS"
            if checks[
                "valid_repeat_counts_match"
            ]
            else "FAIL"
        )
    )

    print(
        "Seed schedule: "
        + (
            "PASS"
            if checks[
                "seed_schedule_match"
            ]
            else "FAIL"
        )
    )

    print(
        "LOO keys: "
        + (
            "PASS"
            if checks[
                "loo_keys_match"
            ]
            else "FAIL"
        )
    )

    print(
        "LOO correlations: "
        + (
            "PASS"
            if checks[
                "loo_correlations_match"
            ]
            else "FAIL"
        )
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
        "Docker backend: "
        f"{backend.kind}"
    )

    print(
        "Image ID: "
        f"{current_image_id}"
    )

    print(
        "Manifest: "
        + str(
            OUTPUT_MANIFEST
            .relative_to(ROOT)
        )
    )

    print("=" * 80)


if __name__ == "__main__":
    main()
