from __future__ import annotations

import ast
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


TOOL_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]

SOURCE_ROOT = (
    ROOT.parent
    / "paper4_sarcoma_dog"
)

CONFIG_PATH = (
    ROOT
    / "configs"
    / "osteosarcoma_reference_registry.json"
)

OUTPUT_ROOT = (
    ROOT
    / "reference_results"
    / "osteosarcoma_locked"
)

ARTIFACT_DIR = (
    OUTPUT_ROOT
    / "artifacts"
)

MANIFEST_DIR = (
    OUTPUT_ROOT
    / "manifests"
)

REGISTRY_PATH = (
    OUTPUT_ROOT
    / "registry.json"
)

SUMMARY_PATH = (
    OUTPUT_ROOT
    / "README.md"
)


VERSION_NAMES = {
    "SCRIPT_VERSION",
    "VERSION",
    "ANALYSIS_VERSION",
}


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


def run_git(
    args: list[str],
) -> str:
    completed = subprocess.run(
        ["git", *args],
        cwd=SOURCE_ROOT,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=True,
    )

    return completed.stdout.strip()


def is_git_tracked(
    relative_path: str,
) -> bool:
    completed = subprocess.run(
        [
            "git",
            "ls-files",
            "--error-unmatch",
            relative_path,
        ],
        cwd=SOURCE_ROOT,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )

    return completed.returncode == 0


def path_git_status(
    relative_path: str,
) -> str:
    return run_git(
        [
            "status",
            "--porcelain=v1",
            "--",
            relative_path,
        ]
    )


def path_last_commit(
    relative_path: str,
) -> dict[str, str]:
    output = run_git(
        [
            "log",
            "-1",
            "--format=%H%x09%cI%x09%s",
            "--",
            relative_path,
        ]
    )

    if not output:
        return {
            "commit": "",
            "date": "",
            "subject": "",
        }

    parts = output.split(
        "\t",
        2,
    )

    while len(parts) < 3:
        parts.append("")

    return {
        "commit": parts[0],
        "date": parts[1],
        "subject": parts[2],
    }


def extract_script_version(
    path: Path,
) -> str:
    text = path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    try:
        tree = ast.parse(text)
    except SyntaxError:
        return ""

    for node in tree.body:
        if not isinstance(
            node,
            (
                ast.Assign,
                ast.AnnAssign,
            ),
        ):
            continue

        if isinstance(
            node,
            ast.Assign,
        ):
            targets = node.targets
            value = node.value

        else:
            targets = [node.target]
            value = node.value

        for target in targets:
            if (
                isinstance(
                    target,
                    ast.Name,
                )
                and target.id
                in VERSION_NAMES
                and isinstance(
                    value,
                    ast.Constant,
                )
                and isinstance(
                    value.value,
                    str,
                )
            ):
                return value.value

    return ""


def require_clean_tracked_path(
    relative_path: str,
) -> None:
    if not is_git_tracked(
        relative_path
    ):
        raise RuntimeError(
            "Reference source is not tracked by Git:\n"
            f"  {relative_path}\n\n"
            "Commit the authoritative source file "
            "before building the registry."
        )

    status = path_git_status(
        relative_path
    )

    if status:
        raise RuntimeError(
            "Reference source has uncommitted changes:\n"
            f"  {relative_path}\n"
            f"Git status: {status}\n\n"
            "Reference registry requires the selected "
            "source path to match a committed Git state."
        )


def load_config() -> dict:
    return json.loads(
        CONFIG_PATH.read_text(
            encoding="utf-8"
        )
    )


def prepare_output_dirs() -> None:
    ARTIFACT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    MANIFEST_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )


def copy_reference_artifact(
    artifact_id: str,
    source_path: Path,
) -> Path:
    destination = (
        ARTIFACT_DIR
        / f"{artifact_id}__{source_path.name}"
    )

    shutil.copy2(
        source_path,
        destination,
    )

    source_sha = sha256_file(
        source_path
    )

    copied_sha = sha256_file(
        destination
    )

    if source_sha != copied_sha:
        raise RuntimeError(
            "Reference copy hash mismatch:\n"
            f"source={source_sha}\n"
            f"copy={copied_sha}"
        )

    return destination


def build_artifact_record(
    definition: dict,
    source_head: str,
) -> dict:
    artifact_id = definition["id"]

    artifact_relative = definition[
        "source_path"
    ]

    script_relative = definition[
        "source_script"
    ]

    artifact_path = (
        SOURCE_ROOT
        / artifact_relative
    )

    script_path = (
        SOURCE_ROOT
        / script_relative
    )

    if not artifact_path.exists():
        raise FileNotFoundError(
            f"Missing artifact: "
            f"{artifact_path}"
        )

    if not script_path.exists():
        raise FileNotFoundError(
            f"Missing source script: "
            f"{script_path}"
        )

    require_clean_tracked_path(
        artifact_relative
    )

    require_clean_tracked_path(
        script_relative
    )

    observed_version = (
        extract_script_version(
            script_path
        )
    )

    expected_version = (
        definition.get(
            "expected_script_version",
            "",
        )
    )

    if (
        expected_version
        and observed_version
        != expected_version
    ):
        raise RuntimeError(
            "Script version mismatch:\n"
            f"script={script_relative}\n"
            f"expected={expected_version}\n"
            f"observed={observed_version}"
        )

    artifact_sha = sha256_file(
        artifact_path
    )

    script_sha = sha256_file(
        script_path
    )

    artifact_git = path_last_commit(
        artifact_relative
    )

    script_git = path_last_commit(
        script_relative
    )

    destination = copy_reference_artifact(
        artifact_id,
        artifact_path,
    )

    record = {
        "id": artifact_id,
        "axis": definition["axis"],
        "regression_class": (
            definition[
                "regression_class"
            ]
        ),
        "source_repository": (
            SOURCE_ROOT.name
        ),
        "source_repository_head_at_registry_build": (
            source_head
        ),
        "source_artifact": {
            "path": artifact_relative,
            "sha256": artifact_sha,
            "size_bytes": (
                artifact_path.stat().st_size
            ),
            "git_commit": (
                artifact_git["commit"]
            ),
            "git_commit_date": (
                artifact_git["date"]
            ),
            "git_commit_subject": (
                artifact_git["subject"]
            ),
        },
        "source_script": {
            "path": script_relative,
            "declared_version": (
                observed_version
            ),
            "sha256": script_sha,
            "git_commit": (
                script_git["commit"]
            ),
            "git_commit_date": (
                script_git["date"]
            ),
            "git_commit_subject": (
                script_git["subject"]
            ),
        },
        "registry_copy": {
            "path": (
                destination.relative_to(
                    ROOT
                ).as_posix()
            ),
            "sha256": sha256_file(
                destination
            ),
            "size_bytes": (
                destination.stat().st_size
            ),
        },
    }

    manifest_path = (
        MANIFEST_DIR
        / f"{artifact_id}.json"
    )

    manifest_path.write_text(
        json.dumps(
            record,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    return record


def detect_duplicate_hashes(
    records: list[dict],
) -> dict[str, list[str]]:
    by_hash: dict[
        str,
        list[str],
    ] = {}

    for record in records:
        sha = record[
            "source_artifact"
        ]["sha256"]

        by_hash.setdefault(
            sha,
            [],
        ).append(
            record["id"]
        )

    return {
        sha: ids
        for sha, ids in by_hash.items()
        if len(ids) > 1
    }


def write_summary(
    registry: dict,
) -> None:
    lines = [
        "# Osteosarcoma Locked Reference Registry",
        "",
        (
            "This directory contains exact copies "
            "of selected authoritative outputs from "
            "the frozen canine-human osteosarcoma study."
        ),
        "",
        (
            "These artifacts are regression targets "
            "for the reusable Molecular Transport Audit "
            "implementation."
        ),
        "",
        f"- Registry version: `{registry['registry_version']}`",
        (
            "- Source repository HEAD at build: "
            f"`{registry['source_repository_head']}`"
        ),
        (
            f"- Artifact count: "
            f"**{len(registry['artifacts'])}**"
        ),
        "",
        "## Artifacts",
        "",
        "| ID | Axis | Regression class | Source SHA-256 |",
        "|---|---|---|---|",
    ]

    for record in registry[
        "artifacts"
    ]:
        lines.append(
            "| "
            f"`{record['id']}` | "
            f"{record['axis']} | "
            f"{record['regression_class']} | "
            f"`{record['source_artifact']['sha256'][:12]}` |"
        )

    duplicates = registry[
        "duplicate_content_hashes"
    ]

    lines.extend(
        [
            "",
            "## Duplicate-content audit",
            "",
        ]
    )

    if not duplicates:
        lines.append(
            "No duplicate reference artifacts detected."
        )

    else:
        for sha, ids in duplicates.items():
            lines.append(
                f"- `{sha}`: "
                + ", ".join(
                    f"`{item}`"
                    for item in ids
                )
            )

    lines.extend(
        [
            "",
            "## Important",
            "",
            (
                "A reference artifact is an immutable "
                "test oracle. It is not automatically "
                "a reusable software input."
            ),
            "",
        ]
    )

    SUMMARY_PATH.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def main() -> None:
    if not SOURCE_ROOT.exists():
        raise FileNotFoundError(
            f"Source repository missing: "
            f"{SOURCE_ROOT}"
        )

    prepare_output_dirs()

    config = load_config()

    source_head = run_git(
        ["rev-parse", "HEAD"]
    )

    records = []

    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- build locked reference registry"
    )
    print("=" * 80)

    print(
        f"Source HEAD: "
        f"{source_head}"
    )

    print(
        f"Definitions: "
        f"{len(config['artifacts'])}"
    )

    print()

    for index, definition in enumerate(
        config["artifacts"],
        start=1,
    ):
        print(
            f"[{index:02d}/"
            f"{len(config['artifacts']):02d}] "
            f"{definition['id']}"
        )

        record = build_artifact_record(
            definition,
            source_head,
        )

        records.append(record)

    duplicates = detect_duplicate_hashes(
        records
    )

    registry = {
        "registry_builder_version": (
            TOOL_VERSION
        ),
        "registry_version": config[
            "registry_version"
        ],
        "study": config["study"],
        "generated_at_utc": (
            datetime.now(
                timezone.utc
            ).isoformat()
        ),
        "source_repository": str(
            SOURCE_ROOT
        ),
        "source_repository_head": (
            source_head
        ),
        "artifacts": records,
        "duplicate_content_hashes": (
            duplicates
        ),
    }

    REGISTRY_PATH.write_text(
        json.dumps(
            registry,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    write_summary(registry)

    print()
    print(
        f"Registry written: "
        f"{REGISTRY_PATH.relative_to(ROOT)}"
    )

    print(
        f"Reference artifacts: "
        f"{len(records)}"
    )

    print(
        f"Duplicate hashes: "
        f"{len(duplicates)}"
    )

    print("=" * 80)


if __name__ == "__main__":
    main()
