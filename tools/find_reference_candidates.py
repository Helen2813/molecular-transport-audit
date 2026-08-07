from __future__ import annotations

import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path


TOOL_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = ROOT.parent / "paper4_sarcoma_dog"

RESULTS_TABLES = SOURCE_ROOT / "results" / "tables"
DOCS_DIR = ROOT / "docs"

OUTPUT_CSV = DOCS_DIR / "reference_candidate_inventory.csv"
OUTPUT_MD = DOCS_DIR / "reference_candidate_inventory.md"
OUTPUT_JSON = DOCS_DIR / "reference_candidate_inventory.json"


AUTHORITATIVE_HINTS = [
    "locked",
    "final",
    "manifest",
    "primary",
    "triangulation",
    "preservation",
    "multiplicity",
    "six_dog",
    "six-dog",
    "rediscovery",
]

EXCLUDE_HINTS = [
    "sentence",
    "caption",
    "outline",
    "readme",
]

ALLOWED_SUFFIXES = {
    ".csv",
    ".json",
    ".txt",
    ".tsv",
    ".tex",
}

SCRIPT_VERSION_PATTERNS = [
    re.compile(
        r'"(?:script_version|analysis_version|version)"\s*:\s*"([^"]+)"',
        re.IGNORECASE,
    ),
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def classify_candidate(path: Path) -> str:
    name = path.name.lower()

    positive = sum(
        1
        for hint in AUTHORITATIVE_HINTS
        if hint in name
    )

    negative = sum(
        1
        for hint in EXCLUDE_HINTS
        if hint in name
    )

    if "manifest" in name:
        return "high"

    if positive >= 2 and negative == 0:
        return "high"

    if positive >= 1 and negative == 0:
        return "medium"

    return "low"


def detect_category(path: Path) -> str:
    name = path.name.lower()

    if "gse239948" in name:
        if "wgcna" in name:
            return "external_canine_wgcna"

        if "rediscovery" in name or "blind" in name:
            return "external_canine_blind_rediscovery"

        return "external_canine_direct"

    if "ammons" in name or "single_cell" in name:
        return "single_cell"

    if "necrosis" in name:
        return "necrosis_exploratory"

    if "multiplicity" in name:
        return "multiplicity"

    if "human" in name or "primary_outcomes" in name:
        return "human_outcome_transport"

    if "module_interpretation" in name:
        return "final_evidence_lock"

    if "multidimensional_transport" in name:
        return "final_evidence_lock"

    if "analysis_lock" in name:
        return "final_evidence_lock"

    if "triangulation" in name:
        return "external_canine_triangulation"

    return "other"


def extract_json_metadata(
    path: Path,
) -> dict[str, str]:
    if path.suffix.lower() != ".json":
        return {
            "declared_script_version": "",
            "declared_script": "",
        }

    try:
        payload = json.loads(
            path.read_text(
                encoding="utf-8",
                errors="replace",
            )
        )
    except (json.JSONDecodeError, OSError):
        return {
            "declared_script_version": "",
            "declared_script": "",
        }

    if not isinstance(payload, dict):
        return {
            "declared_script_version": "",
            "declared_script": "",
        }

    version_keys = [
        "script_version",
        "analysis_version",
        "version",
        "pipeline_version",
    ]

    script_keys = [
        "script",
        "source_script",
        "script_name",
        "analysis_script",
    ]

    version = ""

    for key in version_keys:
        value = payload.get(key)

        if value is not None:
            version = str(value)
            break

    script = ""

    for key in script_keys:
        value = payload.get(key)

        if value is not None:
            script = str(value)
            break

    return {
        "declared_script_version": version,
        "declared_script": script,
    }


def inventory_row(path: Path) -> dict[str, object]:
    relative = path.relative_to(
        SOURCE_ROOT
    ).as_posix()

    metadata = extract_json_metadata(path)

    return {
        "relative_path": relative,
        "filename": path.name,
        "suffix": path.suffix.lower(),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "category": detect_category(path),
        "candidate_priority": classify_candidate(path),
        "declared_script": metadata["declared_script"],
        "declared_script_version": (
            metadata["declared_script_version"]
        ),
        "selected_for_reference_registry": "",
        "notes": "",
    }


def discover_files() -> list[Path]:
    if not RESULTS_TABLES.exists():
        raise FileNotFoundError(
            f"Missing source directory: {RESULTS_TABLES}"
        )

    files = [
        path
        for path in RESULTS_TABLES.iterdir()
        if (
            path.is_file()
            and path.suffix.lower()
            in ALLOWED_SUFFIXES
        )
    ]

    return sorted(
        files,
        key=lambda path: path.name.lower(),
    )


def write_csv(
    rows: list[dict[str, object]],
) -> None:
    fieldnames = list(rows[0].keys())

    with OUTPUT_CSV.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(rows)


def write_markdown(
    rows: list[dict[str, object]],
) -> None:
    high = [
        row
        for row in rows
        if row["candidate_priority"] == "high"
    ]

    medium = [
        row
        for row in rows
        if row["candidate_priority"] == "medium"
    ]

    lines = [
        "# Reference Candidate Inventory",
        "",
        (
            "> Automated candidate discovery only. "
            "No result becomes authoritative merely "
            "because it appears in this report."
        ),
        "",
        f"- Generated: `{datetime.now(timezone.utc).isoformat()}`",
        f"- Total candidate files scanned: **{len(rows)}**",
        f"- High-priority candidates: **{len(high)}**",
        f"- Medium-priority candidates: **{len(medium)}**",
        "",
        "## High-priority candidates",
        "",
        "| Category | File | Size | Script version | SHA-256 |",
        "|---|---|---:|---|---|",
    ]

    for row in high:
        lines.append(
            "| "
            f"{row['category']} | "
            f"`{row['filename']}` | "
            f"{row['size_bytes']} | "
            f"`{row['declared_script_version']}` | "
            f"`{str(row['sha256'])[:12]}` |"
        )

    lines.extend(
        [
            "",
            "## Medium-priority candidates",
            "",
            "| Category | File | Size | Script version |",
            "|---|---|---:|---|",
        ]
    )

    for row in medium:
        lines.append(
            "| "
            f"{row['category']} | "
            f"`{row['filename']}` | "
            f"{row['size_bytes']} | "
            f"`{row['declared_script_version']}` |"
        )

    lines.extend(
        [
            "",
            "## Rule",
            "",
            (
                "Only manually approved files will be copied "
                "into the standalone reference registry."
            ),
            "",
        ]
    )

    OUTPUT_MD.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def write_json(
    rows: list[dict[str, object]],
) -> None:
    payload = {
        "tool_version": TOOL_VERSION,
        "generated_at_utc": (
            datetime.now(
                timezone.utc
            ).isoformat()
        ),
        "source_repository": str(
            SOURCE_ROOT
        ),
        "candidate_count": len(rows),
        "candidates": rows,
    }

    OUTPUT_JSON.write_text(
        json.dumps(
            payload,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    DOCS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    files = discover_files()

    if not files:
        raise RuntimeError(
            "No candidate result files found."
        )

    rows = [
        inventory_row(path)
        for path in files
    ]

    write_csv(rows)
    write_markdown(rows)
    write_json(rows)

    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- reference candidate discovery"
    )
    print("=" * 80)

    print(
        f"Source tables: {RESULTS_TABLES}"
    )

    print(
        f"Files scanned: {len(rows)}"
    )

    high = sum(
        row["candidate_priority"] == "high"
        for row in rows
    )

    medium = sum(
        row["candidate_priority"] == "medium"
        for row in rows
    )

    print(
        f"High priority: {high}"
    )

    print(
        f"Medium priority: {medium}"
    )

    print()

    for row in rows:
        if row["candidate_priority"] != "high":
            continue

        print(
            f"[HIGH] "
            f"{row['category']:35s} "
            f"{row['filename']}"
        )

    print()
    print("Written:")
    print(
        f"  {OUTPUT_CSV.relative_to(ROOT)}"
    )
    print(
        f"  {OUTPUT_MD.relative_to(ROOT)}"
    )
    print(
        f"  {OUTPUT_JSON.relative_to(ROOT)}"
    )

    print("=" * 80)


if __name__ == "__main__":
    main()
