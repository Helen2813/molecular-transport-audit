"""Validate direct preservation inside the Docker scientific environment."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


RUNNER_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]

IMAGE_TAG = (
    "molecular-transport-audit-python:0.1.0"
)

DOCKERFILE = (
    ROOT
    / "containers"
    / "python-core"
    / "Dockerfile"
)

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
    / "docker_preservation_validation"
)

VERIFICATION_FILE = (
    REPORT_ROOT
    / "preservation_regression"
    / "GSE239948_preservation_verification.json"
)

OUTPUT_MANIFEST = (
    REPORT_ROOT
    / "docker_preservation_validation_manifest.json"
)


@dataclass(frozen=True)
class DockerBackend:
    """Docker execution backend."""

    kind: str
    launcher: tuple[str, ...]
    description: str
    wsl_exe: str | None = None

    def command(
        self,
        args: list[str],
    ) -> list[str]:
        return [
            *self.launcher,
            *args,
        ]


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


def run_host(
    command: list[str],
    *,
    capture: bool = False,
    check: bool = True,
    timeout: int | None = None,
) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=(
            subprocess.PIPE
            if capture
            else None
        ),
        stderr=(
            subprocess.PIPE
            if capture
            else None
        ),
        check=False,
        timeout=timeout,
    )

    if (
        check
        and completed.returncode != 0
    ):
        if capture:
            if completed.stdout:
                print(
                    completed.stdout
                )

            if completed.stderr:
                print(
                    completed.stderr
                )

        raise RuntimeError(
            "Command failed with "
            f"exit code "
            f"{completed.returncode}:\n"
            + " ".join(command)
        )

    return completed


def print_command(
    command: list[str],
) -> None:
    printable = " ".join(
        (
            f'"{part}"'
            if " " in str(part)
            else str(part)
        )
        for part in command
    )

    print(
        f"$ {printable}"
    )


def run_docker(
    backend: DockerBackend,
    args: list[str],
    *,
    capture: bool = False,
    check: bool = True,
    timeout: int | None = None,
) -> subprocess.CompletedProcess[str]:
    command = backend.command(
        args
    )

    print_command(
        command
    )

    return run_host(
        command,
        capture=capture,
        check=check,
        timeout=timeout,
    )


def native_docker_candidates(
) -> list[Path]:
    candidates: list[Path] = []

    for executable in [
        "docker",
        "docker.exe",
    ]:
        detected = shutil.which(
            executable
        )

        if detected:
            candidates.append(
                Path(detected)
            )

    for variable in [
        "ProgramFiles",
        "ProgramW6432",
    ]:
        value = os.environ.get(
            variable
        )

        if not value:
            continue

        candidates.append(
            Path(value)
            / "Docker"
            / "Docker"
            / "resources"
            / "bin"
            / "docker.exe"
        )

    local_app_data = os.environ.get(
        "LOCALAPPDATA"
    )

    if local_app_data:
        candidates.append(
            Path(local_app_data)
            / "Programs"
            / "Docker"
            / "Docker"
            / "resources"
            / "bin"
            / "docker.exe"
        )

    unique: list[Path] = []
    seen: set[str] = set()

    for candidate in candidates:
        key = str(
            candidate
        ).lower()

        if key in seen:
            continue

        seen.add(key)
        unique.append(
            candidate
        )

    return unique


def find_native_backend(
) -> DockerBackend | None:
    for candidate in (
        native_docker_candidates()
    ):
        if not candidate.is_file():
            continue

        completed = run_host(
            [
                str(candidate),
                "--version",
            ],
            capture=True,
            check=False,
            timeout=15,
        )

        if completed.returncode == 0:
            return DockerBackend(
                kind=(
                    "native_windows"
                ),
                launcher=(
                    str(candidate),
                ),
                description=(
                    "Native Windows "
                    "Docker CLI"
                ),
            )

    return None


def find_wsl_executable(
) -> str | None:
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

    return None


def find_wsl_backend(
) -> DockerBackend | None:
    wsl_exe = (
        find_wsl_executable()
    )

    if wsl_exe is None:
        return None

    completed = run_host(
        [
            wsl_exe,
            "-e",
            "bash",
            "-lc",
            "command -v docker",
        ],
        capture=True,
        check=False,
        timeout=20,
    )

    if completed.returncode != 0:
        return None

    docker_path = (
        completed.stdout.strip()
    )

    if not docker_path:
        return None

    check = run_host(
        [
            wsl_exe,
            "-e",
            docker_path,
            "--version",
        ],
        capture=True,
        check=False,
        timeout=20,
    )

    if check.returncode != 0:
        return None

    return DockerBackend(
        kind="wsl",
        launcher=(
            wsl_exe,
            "-e",
            docker_path,
        ),
        description=(
            "Docker CLI inside "
            "default WSL distribution "
            f"({docker_path})"
        ),
        wsl_exe=wsl_exe,
    )


def engine_ready(
    backend: DockerBackend,
) -> tuple[
    bool,
    str,
]:
    completed = run_host(
        backend.command(
            [
                "info",
                "--format",
                "{{.ServerVersion}}",
            ]
        ),
        capture=True,
        check=False,
        timeout=20,
    )

    if completed.returncode != 0:
        return (
            False,
            "",
        )

    version = (
        completed.stdout.strip()
    )

    return (
        bool(version),
        version,
    )


def resolve_backend(
) -> DockerBackend:
    print("")
    print(
        "Detecting Docker "
        "execution backend..."
    )

    native = (
        find_native_backend()
    )

    if native is not None:
        ready, version = (
            engine_ready(
                native
            )
        )

        if ready:
            print(
                "Native Docker engine "
                f"ready: {version}"
            )

            return native

    wsl_backend = (
        find_wsl_backend()
    )

    if wsl_backend is not None:
        print(
            "WSL Docker CLI found:"
        )
        print(
            "  "
            + wsl_backend.description
        )

        ready, version = (
            engine_ready(
                wsl_backend
            )
        )

        if ready:
            print(
                "WSL Docker engine "
                f"ready: {version}"
            )

            return wsl_backend

    raise RuntimeError(
        "No usable Docker engine "
        "was found.\n"
        "The preservation implementation "
        "was not executed or modified."
    )


def backend_path(
    backend: DockerBackend,
    path: Path,
) -> str:
    resolved = path.resolve()

    if (
        backend.kind
        == "native_windows"
    ):
        return str(resolved)

    if (
        backend.kind != "wsl"
        or backend.wsl_exe
        is None
    ):
        raise RuntimeError(
            "Cannot convert path for "
            f"backend {backend.kind}."
        )

    completed = run_host(
        [
            backend.wsl_exe,
            "-e",
            "wslpath",
            "-a",
            "-u",
            str(resolved),
        ],
        capture=True,
        timeout=20,
    )

    converted = (
        completed.stdout.strip()
    )

    if not converted:
        raise RuntimeError(
            "WSL path conversion "
            "returned an empty path."
        )

    return converted


def bind_mount(
    backend: DockerBackend,
    source: Path,
    target: str,
    *,
    readonly: bool,
) -> str:
    source_text = backend_path(
        backend,
        source,
    )

    fields = [
        "type=bind",
        f"source={source_text}",
        f"target={target}",
    ]

    if readonly:
        fields.append(
            "readonly"
        )

    return ",".join(
        fields
    )


def docker_versions(
    backend: DockerBackend,
) -> dict[str, str]:
    completed = run_docker(
        backend,
        [
            "version",
            "--format",
            (
                "{{.Client.Version}}"
                "|"
                "{{.Server.Version}}"
            ),
        ],
        capture=True,
    )

    raw = (
        completed.stdout.strip()
    )

    if "|" not in raw:
        return {
            "client": "",
            "server": "",
            "raw": raw,
        }

    client, server = (
        raw.split(
            "|",
            1,
        )
    )

    return {
        "client": client,
        "server": server,
        "raw": raw,
    }


def build_image(
    backend: DockerBackend,
) -> None:
    print("")
    print("=" * 80)
    print(
        "Build preservation-enabled "
        "Python scientific image"
    )
    print("=" * 80)

    dockerfile = backend_path(
        backend,
        DOCKERFILE,
    )

    context = backend_path(
        backend,
        ROOT,
    )

    run_docker(
        backend,
        [
            "build",
            "-f",
            dockerfile,
            "-t",
            IMAGE_TAG,
            context,
        ],
    )


def image_id(
    backend: DockerBackend,
) -> str:
    completed = run_docker(
        backend,
        [
            "image",
            "inspect",
            "--format",
            "{{.Id}}",
            IMAGE_TAG,
        ],
        capture=True,
    )

    return (
        completed.stdout.strip()
    )


def container_versions(
    backend: DockerBackend,
) -> dict[str, str]:
    code = (
        "import json,sys;"
        "import numpy,pandas,yaml,scipy,sklearn;"
        "import transport_audit;"
        "print(json.dumps({"
        "'python':sys.version.split()[0],"
        "'numpy':numpy.__version__,"
        "'pandas':pandas.__version__,"
        "'scipy':scipy.__version__,"
        "'scikit_learn':sklearn.__version__,"
        "'pyyaml':yaml.__version__,"
        "'transport_audit':"
        "transport_audit.__version__"
        "}))"
    )

    completed = run_docker(
        backend,
        [
            "run",
            "--rm",
            "--network",
            "none",
            IMAGE_TAG,
            "python3",
            "-c",
            code,
        ],
        capture=True,
    )

    lines = [
        line.strip()
        for line
        in completed.stdout.splitlines()
        if line.strip()
    ]

    if not lines:
        raise RuntimeError(
            "Container version command "
            "returned no output."
        )

    return json.loads(
        lines[-1]
    )


def run_validation(
    backend: DockerBackend,
) -> None:
    print("")
    print("=" * 80)
    print(
        "Run direct preservation "
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
            "run_preservation_validation_suite.py"
        ),
    ]

    run_docker(
        backend,
        command,
    )


def load_verification(
) -> dict:
    if not VERIFICATION_FILE.exists():
        raise FileNotFoundError(
            "Container preservation "
            "verification was not created:\n"
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
            "Container preservation "
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
    comparison = verification[
        "comparison"
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
            (
                "deterministic_direct_"
                "preservation"
            ),
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
            "dockerfile":
                str(
                    DOCKERFILE
                    .relative_to(ROOT)
                ),
            "dockerfile_sha256":
                sha256_file(
                    DOCKERFILE
                ),
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
        "preservation_verification": {
            "path":
                str(
                    VERIFICATION_FILE
                    .relative_to(ROOT)
                ),
            "sha256":
                sha256_file(
                    VERIFICATION_FILE
                ),
            "passed":
                bool(
                    verification[
                        "passed"
                    ]
                ),
            "fixture_hashes_match":
                bool(
                    verification[
                        "fixture_hashes_match"
                    ]
                ),
            "direct_metrics_match":
                bool(
                    comparison[
                        "direct_metrics_match"
                    ]
                ),
            "coverage_match":
                bool(
                    comparison[
                        "coverage_match"
                    ]
                ),
            "max_absolute_metric_difference":
                float(
                    comparison[
                        "max_absolute_metric_difference"
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
        DOCKERFILE,
        "Dockerfile",
    )

    require_path(
        (
            TOOLS_DIR
            / "run_preservation_validation_suite.py"
        ),
        "Preservation validation suite",
    )

    require_path(
        (
            TESTS_DIR
            / "unit"
            / "test_preservation.py"
        ),
        "Preservation unit tests",
    )

    require_path(
        (
            FIXTURE_DIR
            / "fixture_manifest.json"
        ),
        "Preservation fixture manifest",
    )

    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- Docker preservation validation"
    )
    print("=" * 80)
    print(
        f"Runner version: "
        f"{RUNNER_VERSION}"
    )
    print(
        f"Image: {IMAGE_TAG}"
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

    run_validation(
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

    comparison = verification[
        "comparison"
    ]

    dimensions = verification[
        "dimensions"
    ]

    print("")
    print("=" * 80)
    print(
        "Docker preservation "
        "validation: PASS"
    )
    print("=" * 80)

    print(
        "Fixture hashes: PASS"
    )

    print(
        "Unit and edge-case tests: PASS"
    )

    print(
        "Direct preservation "
        "regression: PASS"
    )

    print(
        "Coverage regression: PASS"
    )

    print(
        "Reference samples: "
        f"{dimensions['reference_samples']}"
    )

    print(
        "External samples: "
        f"{dimensions['external_samples']}"
    )

    print(
        "Fixture genes: "
        f"{dimensions['fixture_genes']}"
    )

    print(
        "Maximum absolute metric "
        "difference: "
        f"{comparison['max_absolute_metric_difference']:.3e}"
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
