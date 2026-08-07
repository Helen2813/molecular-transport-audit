from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


TOOL_VERSION = "0.1.0"

NEW_REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_REPO_ROOT = NEW_REPO_ROOT.parent / "paper4_sarcoma_dog"

DOCS_DIR = NEW_REPO_ROOT / "docs"

JSON_OUTPUT = DOCS_DIR / "source_repository_snapshot.json"
MD_OUTPUT = DOCS_DIR / "source_repository_snapshot.md"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def run_git(
    args: list[str],
    check: bool = True,
) -> str:
    completed = subprocess.run(
        ["git", *args],
        cwd=SOURCE_REPO_ROOT,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=check,
    )

    return completed.stdout.rstrip()


def parse_porcelain_status(
    text: str,
) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []

    for raw_line in text.splitlines():
        if len(raw_line) < 4:
            continue

        index_status = raw_line[0]
        worktree_status = raw_line[1]

        path_text = raw_line[3:]

        original_path = ""
        current_path = path_text

        if " -> " in path_text:
            original_path, current_path = (
                path_text.split(" -> ", 1)
            )

        entries.append(
            {
                "index_status": index_status,
                "worktree_status": worktree_status,
                "path": current_path,
                "original_path": original_path,
            }
        )

    return entries


def classify_status(
    entry: dict[str, str],
) -> str:
    index_status = entry["index_status"]
    worktree_status = entry["worktree_status"]

    if index_status == "?" and worktree_status == "?":
        return "untracked"

    statuses = {
        index_status,
        worktree_status,
    }

    if "D" in statuses:
        return "deleted"

    if "R" in statuses:
        return "renamed"

    if "A" in statuses:
        return "added"

    if "M" in statuses:
        return "modified"

    return "other"


def snapshot_path(
    relative_path: str,
) -> dict[str, object]:
    full_path = SOURCE_REPO_ROOT / relative_path

    result: dict[str, object] = {
        "path": relative_path,
        "exists": full_path.exists(),
        "is_file": full_path.is_file(),
        "size_bytes": None,
        "sha256": None,
    }

    if full_path.is_file():
        result["size_bytes"] = full_path.stat().st_size
        result["sha256"] = sha256_file(full_path)

    return result


def source_metadata() -> dict[str, object]:
    head = run_git(
        ["rev-parse", "HEAD"]
    )

    branch = run_git(
        ["branch", "--show-current"]
    )

    remote = run_git(
        ["remote", "get-url", "origin"],
        check=False,
    )

    status_text = run_git(
        ["status", "--porcelain=v1"]
    )

    diff_name_status = run_git(
        ["diff", "--name-status"]
    )

    diff_stat = run_git(
        ["diff", "--stat"]
    )

    staged_diff_name_status = run_git(
        [
            "diff",
            "--cached",
            "--name-status",
        ]
    )

    entries = parse_porcelain_status(
        status_text
    )

    file_snapshots = []

    for entry in entries:
        current_path = entry["path"]

        file_snapshot = snapshot_path(
            current_path
        )

        file_snapshot.update(
            {
                "status_class": (
                    classify_status(entry)
                ),
                "index_status": (
                    entry["index_status"]
                ),
                "worktree_status": (
                    entry["worktree_status"]
                ),
                "original_path": (
                    entry["original_path"]
                ),
            }
        )

        file_snapshots.append(
            file_snapshot
        )

    return {
        "snapshot_tool_version": TOOL_VERSION,
        "generated_at_utc": (
            datetime.now(
                timezone.utc
            ).isoformat()
        ),
        "source_repository": (
            SOURCE_REPO_ROOT.name
        ),
        "source_repository_path": str(
            SOURCE_REPO_ROOT
        ),
        "remote": remote,
        "branch": branch,
        "head_commit": head,
        "worktree_dirty": bool(entries),
        "status_entries": entries,
        "changed_files": file_snapshots,
        "unstaged_diff_name_status": (
            diff_name_status.splitlines()
            if diff_name_status
            else []
        ),
        "staged_diff_name_status": (
            staged_diff_name_status.splitlines()
            if staged_diff_name_status
            else []
        ),
        "diff_stat": (
            diff_stat.splitlines()
            if diff_stat
            else []
        ),
    }


def write_markdown(
    snapshot: dict[str, object],
) -> None:
    changed_files = snapshot[
        "changed_files"
    ]

    lines = [
        "# Source Repository Snapshot",
        "",
        (
            "> This snapshot records the local state "
            "of the frozen osteosarcoma source repository "
            "before reusable software migration."
        ),
        "",
        "## Repository",
        "",
        (
            f"- Repository: "
            f"`{snapshot['source_repository']}`"
        ),
        (
            f"- Path: "
            f"`{snapshot['source_repository_path']}`"
        ),
        (
            f"- Remote: "
            f"`{snapshot['remote']}`"
        ),
        (
            f"- Branch: "
            f"`{snapshot['branch']}`"
        ),
        (
            f"- HEAD commit: "
            f"`{snapshot['head_commit']}`"
        ),
        (
            f"- Worktree dirty: "
            f"**{snapshot['worktree_dirty']}**"
        ),
        (
            f"- Generated: "
            f"`{snapshot['generated_at_utc']}`"
        ),
        "",
        "## Changed files",
        "",
    ]

    if not changed_files:
        lines.append(
            "No changed or untracked files detected."
        )

    else:
        lines.extend(
            [
                "| Status | Path | Size | SHA-256 |",
                "|---|---|---:|---|",
            ]
        )

        for file_info in changed_files:
            sha = file_info.get(
                "sha256"
            )

            short_sha = (
                str(sha)[:12]
                if sha
                else ""
            )

            size = file_info.get(
                "size_bytes"
            )

            lines.append(
                "| "
                f"{file_info['status_class']} | "
                f"`{file_info['path']}` | "
                f"{size if size is not None else ''} | "
                f"`{short_sha}` |"
            )

    lines.extend(
        [
            "",
            "## Diff summary",
            "",
            "```text",
            *snapshot["diff_stat"],
            "```",
            "",
            "## Interpretation",
            "",
            (
                "A dirty worktree does not automatically "
                "invalidate the source study."
            ),
            "",
            (
                "However, authoritative reference outputs "
                "must be linked to the exact source files "
                "and versions that generated them."
            ),
            "",
            (
                "No scientific migration should rely only "
                "on the repository HEAD when locally modified "
                "analysis files are present."
            ),
            "",
        ]
    )

    MD_OUTPUT.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def main() -> None:
    if not SOURCE_REPO_ROOT.exists():
        raise FileNotFoundError(
            "Source repository not found:\n"
            f"{SOURCE_REPO_ROOT}"
        )

    DOCS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    snapshot = source_metadata()

    JSON_OUTPUT.write_text(
        json.dumps(
            snapshot,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    write_markdown(snapshot)

    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- source repository snapshot"
    )
    print("=" * 80)

    print(
        f"Source repository: "
        f"{SOURCE_REPO_ROOT}"
    )

    print(
        f"HEAD: "
        f"{snapshot['head_commit']}"
    )

    print(
        f"Branch: "
        f"{snapshot['branch']}"
    )

    print(
        f"Dirty: "
        f"{snapshot['worktree_dirty']}"
    )

    changed_files = snapshot[
        "changed_files"
    ]

    print(
        f"Changed/untracked files: "
        f"{len(changed_files)}"
    )

    print()

    for item in changed_files:
        print(
            f"  {item['status_class']:10s} "
            f"{item['path']}"
        )

    print()

    print(
        "Written:"
    )

    print(
        "  "
        + str(
            JSON_OUTPUT.relative_to(
                NEW_REPO_ROOT
            )
        )
    )

    print(
        "  "
        + str(
            MD_OUTPUT.relative_to(
                NEW_REPO_ROOT
            )
        )
    )

    print("=" * 80)


if __name__ == "__main__":
    main()
