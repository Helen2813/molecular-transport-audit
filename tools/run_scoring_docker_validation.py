"""Validate frozen-program scoring inside Docker on Windows or WSL."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


RUNNER_VERSION = "0.2.0"

ROOT = Path(__file__).resolve().parents[1]

IMAGE_TAG = "molecular-transport-audit-python:0.1.0"

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
    / "scoring_fixture"
)

REPORT_ROOT = (
    ROOT
    / "reports"
    / "docker_scoring_validation"
)

SCORING_REPORT = (
    REPORT_ROOT
    / "scoring_regression"
    / "TARGET_OS_scoring_verification.json"
)

OUTPUT_MANIFEST = (
    REPORT_ROOT
    / "docker_scoring_validation_manifest.json"
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


def run_host_command(
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

    if check and completed.returncode != 0:
        if capture:
            if completed.stdout:
                print(completed.stdout)

            if completed.stderr:
                print(completed.stderr)

        raise RuntimeError(
            "Command failed with "
            f"exit code {completed.returncode}:\n"
            + " ".join(command)
        )

    return completed


def print_command(
    command: list[str],
) -> None:
    printable = " ".join(
        f'"{part}"'
        if " " in str(part)
        else str(part)
        for part in command
    )

    print(f"$ {printable}")


def run_docker(
    backend: DockerBackend,
    args: list[str],
    *,
    capture: bool = False,
    check: bool = True,
    timeout: int | None = None,
) -> subprocess.CompletedProcess[str]:
    command = backend.command(args)

    print_command(command)

    return run_host_command(
        command,
        capture=capture,
        check=check,
        timeout=timeout,
    )


def native_docker_candidates() -> list[Path]:
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

    program_files_values = [
        os.environ.get("ProgramFiles"),
        os.environ.get("ProgramW6432"),
    ]

    for value in program_files_values:
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
        candidates.extend(
            [
                (
                    Path(local_app_data)
                    / "Programs"
                    / "Docker"
                    / "Docker"
                    / "resources"
                    / "bin"
                    / "docker.exe"
                ),
                (
                    Path(local_app_data)
                    / "Docker"
                    / "resources"
                    / "bin"
                    / "docker.exe"
                ),
            ]
        )

    candidates.append(
        Path(
            r"C:\ProgramData\DockerDesktop"
            r"\version-bin\docker.exe"
        )
    )

    unique: list[Path] = []
    seen: set[str] = set()

    for candidate in candidates:
        key = str(candidate).lower()

        if key in seen:
            continue

        seen.add(key)
        unique.append(candidate)

    return unique


def find_native_backend() -> DockerBackend | None:
    for candidate in native_docker_candidates():
        if not candidate.is_file():
            continue

        completed = run_host_command(
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
                kind="native_windows",
                launcher=(
                    str(candidate),
                ),
                description=(
                    "Native Windows Docker CLI"
                ),
            )

    return None


def find_wsl_executable() -> str | None:
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


def find_wsl_backend() -> DockerBackend | None:
    wsl_exe = find_wsl_executable()

    if wsl_exe is None:
        return None

    completed = run_host_command(
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

    version_check = run_host_command(
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

    if version_check.returncode != 0:
        return None

    return DockerBackend(
        kind="wsl",
        launcher=(
            wsl_exe,
            "-e",
            docker_path,
        ),
        description=(
            "Docker CLI inside default WSL "
            f"distribution ({docker_path})"
        ),
        wsl_exe=wsl_exe,
    )


def docker_engine_ready(
    backend: DockerBackend,
) -> tuple[bool, str]:
    completed = run_host_command(
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
        return False, ""

    server_version = (
        completed.stdout.strip()
    )

    return bool(
        server_version
    ), server_version


def docker_desktop_candidates() -> list[Path]:
    candidates: list[Path] = []

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
            / "Docker Desktop.exe"
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
            / "Docker Desktop.exe"
        )

    return candidates


def try_start_docker_desktop(
    backend: DockerBackend,
    timeout_seconds: int = 120,
) -> bool:
    if backend.kind != "native_windows":
        return False

    desktop_executable: Path | None = None

    for candidate in docker_desktop_candidates():
        if candidate.is_file():
            desktop_executable = candidate
            break

    if desktop_executable is None:
        return False

    print("")
    print(
        "Docker CLI was found but the "
        "Docker engine is not ready."
    )
    print(
        "Starting Docker Desktop:"
    )
    print(
        f"  {desktop_executable}"
    )

    subprocess.Popen(
        [
            str(desktop_executable)
        ],
        cwd=desktop_executable.parent,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    deadline = (
        time.monotonic()
        + timeout_seconds
    )

    while time.monotonic() < deadline:
        ready, version = (
            docker_engine_ready(
                backend
            )
        )

        if ready:
            print(
                "Docker Desktop engine "
                f"is ready: {version}"
            )
            return True

        print(
            "Waiting for Docker engine..."
        )

        time.sleep(3)

    return False


def resolve_backend() -> DockerBackend:
    print("")
    print(
        "Detecting Docker execution "
        "backend..."
    )

    native = find_native_backend()

    if native is not None:
        print(
            "Native Docker CLI found:"
        )
        print(
            f"  {native.launcher[0]}"
        )

        ready, server = (
            docker_engine_ready(
                native
            )
        )

        if ready:
            print(
                "Native Docker engine "
                f"ready: {server}"
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
            f"  {wsl_backend.description}"
        )

        ready, server = (
            docker_engine_ready(
                wsl_backend
            )
        )

        if ready:
            print(
                "WSL Docker engine "
                f"ready: {server}"
            )
            return wsl_backend

    if native is not None:
        if try_start_docker_desktop(
            native
        ):
            return native

    diagnostics: list[str] = []

    if native is None:
        diagnostics.append(
            "- No usable native Windows "
            "docker.exe was found."
        )
    else:
        diagnostics.append(
            "- Native Docker CLI exists, "
            "but its engine is unavailable."
        )

    if wsl_backend is None:
        diagnostics.append(
            "- Docker CLI was not detected "
            "inside the default WSL "
            "distribution."
        )
    else:
        diagnostics.append(
            "- Docker CLI exists inside WSL, "
            "but its Docker engine is "
            "unavailable."
        )

    raise RuntimeError(
        "No usable Docker engine was found.\n\n"
        + "\n".join(diagnostics)
        + "\n\n"
        "The scoring implementation itself "
        "has not been executed or modified."
    )


def backend_path(
    backend: DockerBackend,
    path: Path,
) -> str:
    resolved = path.resolve()

    if backend.kind == "native_windows":
        return str(resolved)

    if backend.kind != "wsl":
        raise RuntimeError(
            f"Unknown backend kind: "
            f"{backend.kind}"
        )

    if not backend.wsl_exe:
        raise RuntimeError(
            "WSL executable was not "
            "recorded."
        )

    completed = run_host_command(
        [
            backend.wsl_exe,
            "-e",
            "wslpath",
            "-a",
            "-u",
            str(resolved),
        ],
        capture=True,
        check=True,
        timeout=20,
    )

    converted = (
        completed.stdout.strip()
    )

    if not converted:
        raise RuntimeError(
            "WSL path conversion returned "
            f"no path for {resolved}"
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

    return ",".join(fields)


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

    text = (
        completed.stdout.strip()
    )

    if "|" not in text:
        return {
            "client": "",
            "server": "",
            "raw": text,
        }

    client, server = (
        text.split("|", 1)
    )

    return {
        "client": client,
        "server": server,
        "raw": text,
    }


def build_image(
    backend: DockerBackend,
) -> None:
    print("")
    print("=" * 80)
    print(
        "Build scoring-enabled "
        "Python scientific image"
    )
    print("=" * 80)

    dockerfile = backend_path(
        backend,
        DOCKERFILE,
    )

    build_context = backend_path(
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
            build_context,
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
        "import json,sys,numpy,pandas,yaml,"
        "transport_audit;"
        "print(json.dumps({"
        "'python':sys.version.split()[0],"
        "'numpy':numpy.__version__,"
        "'pandas':pandas.__version__,"
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


def run_container_validation(
    backend: DockerBackend,
) -> None:
    print("")
    print("=" * 80)
    print(
        "Run frozen scoring suite "
        "inside Docker"
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
                "scoring_fixture"
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
            "run_scoring_validation_suite.py"
        ),
    ]

    run_docker(
        backend,
        command,
    )


def load_verification() -> dict:
    if not SCORING_REPORT.exists():
        raise FileNotFoundError(
            "Container scoring verification "
            "was not created:\n"
            f"{SCORING_REPORT}"
        )

    verification = json.loads(
        SCORING_REPORT.read_text(
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
            "Container scoring verification "
            "reported passed=false."
        )

    return verification


def write_manifest(
    *,
    backend: DockerBackend,
    docker_info: dict[str, str],
    image_sha: str,
    package_versions: dict[str, str],
    verification: dict,
) -> None:
    comparison = verification.get(
        "comparison",
        {},
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
        "scoring_verification": {
            "path":
                str(
                    SCORING_REPORT
                    .relative_to(ROOT)
                ),
            "sha256":
                sha256_file(
                    SCORING_REPORT
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
            "scores_match":
                bool(
                    comparison[
                        "scores_match"
                    ]
                ),
            "coverage_match":
                bool(
                    comparison[
                        "coverage_match"
                    ]
                ),
            "n_samples":
                int(
                    comparison[
                        "n_samples"
                    ]
                ),
            "n_score_columns":
                int(
                    comparison[
                        "n_score_columns"
                    ]
                ),
            "max_absolute_score_difference":
                float(
                    comparison[
                        "max_absolute_score_difference"
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
        TOOLS_DIR
        / "run_scoring_validation_suite.py",
        "Scoring validation suite",
    )

    require_path(
        TESTS_DIR
        / "unit"
        / "test_scoring.py",
        "Scoring unit tests",
    )

    require_path(
        FIXTURE_DIR
        / "fixture_manifest.json",
        "Scoring fixture manifest",
    )

    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- Docker scoring validation"
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

    docker_info = docker_versions(
        backend
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

    current_image_id = image_id(
        backend
    )

    versions = container_versions(
        backend
    )

    print("")
    print(
        "Container environment:"
    )

    for key, value in versions.items():
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

    comparison = verification[
        "comparison"
    ]

    print("")
    print("=" * 80)
    print(
        "Docker scoring validation: PASS"
    )
    print("=" * 80)
    print(
        "Fixture hashes: PASS"
    )
    print(
        "Unit and edge-case tests: PASS"
    )
    print(
        "Score regression: PASS"
    )
    print(
        "Coverage regression: PASS"
    )
    print(
        "Samples: "
        f"{comparison['n_samples']}"
    )
    print(
        "Score columns: "
        f"{comparison['n_score_columns']}"
    )
    print(
        "Maximum absolute score "
        "difference: "
        f"{comparison['max_absolute_score_difference']:.3e}"
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
