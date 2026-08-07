from __future__ import annotations

import ast
import csv
import hashlib
import json
import re
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


TOOL_VERSION = "0.1.0"

# Repositories are expected to be sibling directories.
NEW_REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_REPO_ROOT = NEW_REPO_ROOT.parent / "paper4_sarcoma_dog"
SOURCE_SCRIPTS_DIR = SOURCE_REPO_ROOT / "scripts"

DOCS_DIR = NEW_REPO_ROOT / "docs"

AUTO_CSV = DOCS_DIR / "legacy_analysis_inventory_auto.csv"
AUTO_MD = DOCS_DIR / "legacy_analysis_inventory_auto.md"
REVIEW_CSV = DOCS_DIR / "legacy_analysis_inventory_review.csv"
META_JSON = DOCS_DIR / "legacy_analysis_inventory_meta.json"

SUPPORTED_SUFFIXES = {
    ".py",
    ".r",
    ".rmd",
}

FILE_LITERAL_EXTENSIONS = {
    ".csv",
    ".tsv",
    ".txt",
    ".gz",
    ".xlsx",
    ".xls",
    ".parquet",
    ".feather",
    ".json",
    ".yaml",
    ".yml",
    ".pkl",
    ".pickle",
    ".rds",
    ".h5ad",
    ".h5",
    ".hdf5",
    ".npy",
    ".npz",
    ".png",
    ".jpg",
    ".jpeg",
    ".svg",
    ".pdf",
}

# These fields are reviewed manually and are preserved across rescans.
MANUAL_FIELDS = [
    "authoritative_version",
    "framework_axis",
    "scientific_purpose",
    "reusable_logic",
    "dataset_specific_logic",
    "authoritative_locked_output",
    "proposed_destination_module",
    "regression_strategy",
    "review_status",
    "notes",
]

# These are only automated suggestions.
# They must not be treated as scientific classifications.
AXIS_KEYWORDS = {
    "predictive_specificity": {
        "random control": 4,
        "random panel": 4,
        "random gene": 3,
        "matched random": 4,
        "nested cv": 3,
        "cross-validation": 2,
        "cross validation": 2,
        "c-index": 2,
        "concordance index": 2,
        "elastic net": 2,
        "conditional cox": 2,
        "gene-level": 2,
        "module-level": 2,
    },
    "molecular_representation_preservation": {
        "preservation": 3,
        "edge concordance": 4,
        "loading concordance": 4,
        "split-half": 3,
        "split half": 3,
        "wgcna": 4,
        "zsummary": 4,
        "medianrank": 3,
        "modulepreservation": 4,
        "external canine": 2,
        "network preservation": 3,
    },
    "unsupervised_recurrence": {
        "mofa": 4,
        "mofapy": 4,
        "latent factor": 4,
        "latent-factor": 4,
        "blind rediscovery": 5,
        "de novo": 3,
        "denovo": 3,
        "rediscovery": 4,
        "max-over-factors": 4,
        "factor recurrence": 4,
    },
    "clinical_outcome_transport": {
        "target-os": 3,
        "gse21257": 3,
        "gse39055": 3,
        "outcome transport": 5,
        "overall survival": 3,
        "disease-free interval": 3,
        "recurrence-free survival": 3,
        "hazard ratio": 3,
        "coxph": 2,
        "cox proportional": 2,
        "average precision": 2,
        "metastasis": 2,
    },
    "measurement_robustness": {
        "assay": 3,
        "detection p": 4,
        "detection-aware": 4,
        "detection aware": 4,
        "probe rule": 4,
        "probe-level": 3,
        "ffpe": 3,
        "dasl": 3,
        "measurement rule": 4,
        "quality diagnostic": 3,
    },
    "multiplicity_provenance": {
        "benjamini": 4,
        "hochberg": 4,
        "multiplicity": 4,
        "false discovery rate": 3,
        "q-value": 2,
        "qvalue": 2,
        "manifest": 3,
        "sha256": 3,
        "sha-256": 3,
        "locked": 2,
        "frozen": 1,
        "permutation-adjusted": 3,
    },
}

LOCK_MARKERS = [
    "FIXED_V2",
    "FIXED V2",
    "final lock",
    "locked",
    "frozen",
    "script version",
]

SEED_PATTERNS = [
    re.compile(r"\bseed\s*=\s*(\d+)", re.IGNORECASE),
    re.compile(r"\brandom_state\s*=\s*(\d+)", re.IGNORECASE),
    re.compile(r"default_rng\(\s*(\d+)\s*\)", re.IGNORECASE),
    re.compile(
        r"(?:np|numpy)\.random\.seed\(\s*(\d+)\s*\)",
        re.IGNORECASE,
    ),
    re.compile(r"set\.seed\(\s*(\d+)\s*\)", re.IGNORECASE),
]

VERSION_PATTERNS = [
    re.compile(
        r"^\s*(?:SCRIPT_VERSION|VERSION|ANALYSIS_VERSION)"
        r"\s*=\s*[\"']([^\"']+)[\"']",
        re.MULTILINE,
    ),
    re.compile(
        r"Script version:\s*([A-Za-z0-9._-]+)",
        re.IGNORECASE,
    ),
]

RANDOMIZATION_TERMS = [
    "permutation",
    "bootstrap",
    "resample",
    "random_state",
    "default_rng",
    "np.random",
    "numpy.random",
    "set.seed",
    "shuffle",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def run_git(args: list[str]) -> str:
    try:
        completed = subprocess.run(
            ["git", *args],
            cwd=SOURCE_REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return ""

    return completed.stdout.strip()


def source_repo_metadata() -> dict[str, object]:
    head = run_git(["rev-parse", "HEAD"])
    branch = run_git(["branch", "--show-current"])
    status = run_git(["status", "--porcelain"])
    remote = run_git(["remote", "get-url", "origin"])

    return {
        "inventory_tool_version": TOOL_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_repository_name": SOURCE_REPO_ROOT.name,
        "source_repository_remote": remote,
        "source_head_commit": head,
        "source_branch": branch,
        "source_worktree_dirty": bool(status),
    }


def file_git_metadata(path: Path) -> dict[str, str]:
    relative_path = path.relative_to(SOURCE_REPO_ROOT).as_posix()

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
            "last_git_commit": "",
            "last_git_commit_date": "",
            "last_git_commit_subject": "",
        }

    parts = output.split("\t", 2)

    while len(parts) < 3:
        parts.append("")

    return {
        "last_git_commit": parts[0],
        "last_git_commit_date": parts[1],
        "last_git_commit_subject": parts[2],
    }


def extract_script_number(filename: str) -> int | None:
    match = re.match(r"^(\d+)[_-]", filename)

    if not match:
        return None

    return int(match.group(1))


def sort_key(path: Path) -> tuple[bool, int, str]:
    number = extract_script_number(path.name)

    return (
        number is None,
        number if number is not None else 10**9,
        path.name.lower(),
    )


def read_text(path: Path) -> str:
    return path.read_text(
        encoding="utf-8",
        errors="replace",
    )


def extract_python_details(
    text: str,
) -> tuple[str, list[str], list[str]]:
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return "", [], []

    docstring = ast.get_docstring(tree) or ""

    imports: set[str] = set()
    literals: set[str] = set()

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imports.add(alias.name.split(".")[0])

        elif isinstance(node, ast.ImportFrom):
            if node.module:
                imports.add(node.module.split(".")[0])

        elif isinstance(node, ast.Constant):
            if isinstance(node.value, str):
                value = node.value.strip()

                if value:
                    literals.add(value)

    return (
        docstring,
        sorted(imports),
        sorted(literals),
    )


def extract_r_details(
    text: str,
) -> tuple[str, list[str], list[str]]:
    packages = sorted(
        set(
            re.findall(
                r"(?:library|require)\s*\(\s*[\"']?"
                r"([A-Za-z0-9_.]+)",
                text,
                flags=re.IGNORECASE,
            )
        )
    )

    literals: set[str] = set()

    for match in re.finditer(
        r"([\"'])(.*?)(?<!\\)\1",
        text,
        flags=re.DOTALL,
    ):
        value = match.group(2).strip()

        if value:
            literals.add(value)

    comments: list[str] = []

    for line in text.splitlines():
        stripped = line.strip()

        if stripped.startswith("#"):
            comments.append(
                stripped.lstrip("#").strip()
            )
        elif comments:
            break

    docstring = " ".join(comments[:8])

    return (
        docstring,
        packages,
        sorted(literals),
    )


def looks_like_file_literal(value: str) -> bool:
    cleaned = value.strip().replace("\\", "/")

    if len(cleaned) > 300:
        return False

    lower = cleaned.lower()

    return any(
        lower.endswith(extension)
        for extension in FILE_LITERAL_EXTENSIONS
    )


def extract_file_literals(
    literals: Iterable[str],
) -> list[str]:
    return sorted(
        {
            value
            for value in literals
            if looks_like_file_literal(value)
        }
    )


def extract_version(text: str) -> str:
    for pattern in VERSION_PATTERNS:
        match = pattern.search(text)

        if match:
            return match.group(1).strip()

    if re.search(
        r"\bFIXED_V2\b",
        text,
        flags=re.IGNORECASE,
    ):
        return "FIXED_V2"

    return ""


def extract_seeds(text: str) -> list[str]:
    seeds: set[str] = set()

    for pattern in SEED_PATTERNS:
        seeds.update(pattern.findall(text))

    return sorted(
        seeds,
        key=lambda value: int(value),
    )


def detect_randomization(text: str) -> list[str]:
    lower = text.lower()

    return sorted(
        {
            term
            for term in RANDOMIZATION_TERMS
            if term.lower() in lower
        }
    )


def suggest_axes(
    text: str,
) -> tuple[list[str], dict[str, int]]:
    lower = text.lower()

    scores: dict[str, int] = {}

    for axis, keywords in AXIS_KEYWORDS.items():
        score = 0

        for keyword, weight in keywords.items():
            if keyword.lower() in lower:
                score += weight

        if score > 0:
            scores[axis] = score

    if not scores:
        return [], {}

    max_score = max(scores.values())

    # Require at least moderate evidence before suggesting an axis.
    threshold = max(
        2,
        min(4, max_score),
    )

    selected = [
        axis
        for axis, score in sorted(
            scores.items(),
            key=lambda item: (-item[1], item[0]),
        )
        if score >= threshold
    ]

    return selected, scores


def detect_lock_markers(text: str) -> list[str]:
    lower = text.lower()

    return sorted(
        {
            marker
            for marker in LOCK_MARKERS
            if marker.lower() in lower
        }
    )


def summarize_docstring(
    docstring: str,
    limit: int = 280,
) -> str:
    compact = " ".join(docstring.split())

    if len(compact) <= limit:
        return compact

    return (
        compact[: limit - 3].rstrip()
        + "..."
    )


def scan_script(path: Path) -> dict[str, str]:
    text = read_text(path)

    suffix = path.suffix.lower()

    if suffix == ".py":
        docstring, imports, literals = (
            extract_python_details(text)
        )
        language = "python"

    else:
        docstring, imports, literals = (
            extract_r_details(text)
        )
        language = "r"

    file_literals = extract_file_literals(literals)

    axes, axis_scores = suggest_axes(text)

    git_metadata = file_git_metadata(path)

    number = extract_script_number(path.name)

    line_count = (
        text.count("\n")
        + (1 if text else 0)
    )

    randomization_terms = detect_randomization(text)

    return {
        "script_number": (
            ""
            if number is None
            else str(number)
        ),
        "relative_path": (
            path.relative_to(
                SOURCE_REPO_ROOT
            ).as_posix()
        ),
        "filename": path.name,
        "language": language,
        "declared_version": extract_version(text),
        "sha256": sha256_file(path),
        "size_bytes": str(path.stat().st_size),
        "line_count": str(line_count),
        "last_git_commit": (
            git_metadata["last_git_commit"]
        ),
        "last_git_commit_date": (
            git_metadata["last_git_commit_date"]
        ),
        "last_git_commit_subject": (
            git_metadata[
                "last_git_commit_subject"
            ]
        ),
        "docstring_summary": (
            summarize_docstring(docstring)
        ),
        "imports_or_packages": "; ".join(imports),
        "randomization_detected": (
            "yes"
            if randomization_terms
            else "no"
        ),
        "randomization_terms": (
            "; ".join(randomization_terms)
        ),
        "seed_literals": (
            "; ".join(extract_seeds(text))
        ),
        "candidate_file_literals": (
            "; ".join(file_literals)
        ),
        "candidate_file_literal_count": (
            str(len(file_literals))
        ),
        "auto_axis_suggestion": (
            "; ".join(axes)
        ),
        "auto_axis_scores": json.dumps(
            axis_scores,
            sort_keys=True,
        ),
        "lock_markers": (
            "; ".join(
                detect_lock_markers(text)
            )
        ),
        "same_number_file_count": "",
        "same_number_files": "",
    }


def discover_scripts() -> list[Path]:
    scripts = [
        path
        for path in SOURCE_SCRIPTS_DIR.rglob("*")
        if (
            path.is_file()
            and path.suffix.lower()
            in SUPPORTED_SUFFIXES
        )
    ]

    return sorted(
        scripts,
        key=sort_key,
    )


def annotate_duplicate_numbers(
    rows: list[dict[str, str]],
) -> None:
    groups: dict[str, list[str]] = {}

    for row in rows:
        script_number = row["script_number"]

        if not script_number:
            continue

        groups.setdefault(
            script_number,
            [],
        ).append(row["filename"])

    for row in rows:
        script_number = row["script_number"]

        if not script_number:
            continue

        files = groups.get(
            script_number,
            [],
        )

        row["same_number_file_count"] = str(
            len(files)
        )

        row["same_number_files"] = "; ".join(
            sorted(files)
        )


def write_csv(
    path: Path,
    rows: list[dict[str, str]],
    fieldnames: list[str],
) -> None:
    with path.open(
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


def load_existing_manual_review(
) -> dict[str, dict[str, str]]:
    if not REVIEW_CSV.exists():
        return {}

    existing: dict[
        str,
        dict[str, str],
    ] = {}

    with REVIEW_CSV.open(
        "r",
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        reader = csv.DictReader(handle)

        for row in reader:
            relative_path = row.get(
                "relative_path",
                "",
            )

            if not relative_path:
                continue

            existing[relative_path] = {
                field: row.get(field, "")
                for field in MANUAL_FIELDS
            }

    return existing


def write_review_csv(
    auto_rows: list[dict[str, str]],
    auto_fields: list[str],
) -> None:
    existing = load_existing_manual_review()

    review_rows: list[
        dict[str, str]
    ] = []

    for row in auto_rows:
        merged = dict(row)

        previous_manual_values = existing.get(
            row["relative_path"],
            {},
        )

        for field in MANUAL_FIELDS:
            merged[field] = (
                previous_manual_values.get(
                    field,
                    "",
                )
            )

        review_rows.append(merged)

    write_csv(
        REVIEW_CSV,
        review_rows,
        [
            *auto_fields,
            *MANUAL_FIELDS,
        ],
    )


def short_sha(value: str) -> str:
    if not value:
        return ""

    return value[:10]


def markdown_escape(value: str) -> str:
    return (
        value.replace("|", "\\|")
        .replace("\n", " ")
    )


def write_markdown(
    rows: list[dict[str, str]],
    metadata: dict[str, object],
) -> None:
    axis_counts: Counter[str] = Counter()

    random_count = 0
    version_count = 0
    duplicate_number_count = 0

    seen_duplicate_numbers: set[str] = set()

    for row in rows:
        if row["randomization_detected"] == "yes":
            random_count += 1

        if row["declared_version"]:
            version_count += 1

        if (
            row["same_number_file_count"]
            and int(
                row["same_number_file_count"]
            ) > 1
        ):
            seen_duplicate_numbers.add(
                row["script_number"]
            )

        for axis in filter(
            None,
            (
                part.strip()
                for part in row[
                    "auto_axis_suggestion"
                ].split(";")
            ),
        ):
            axis_counts[axis] += 1

    duplicate_number_count = len(
        seen_duplicate_numbers
    )

    lines = [
        "# Legacy Analysis Inventory — Automated Scan",
        "",
        "> This file is generated automatically. "
        "Axis assignments are suggestions only "
        "and are not scientific classifications.",
        "",
        "## Source snapshot",
        "",
        (
            "- Source repository: "
            f"`{metadata.get('source_repository_name', '')}`"
        ),
        (
            "- Remote: "
            f"`{metadata.get('source_repository_remote', '')}`"
        ),
        (
            "- Branch: "
            f"`{metadata.get('source_branch', '')}`"
        ),
        (
            "- HEAD: "
            f"`{metadata.get('source_head_commit', '')}`"
        ),
        (
            "- Source worktree dirty: "
            f"`{metadata.get('source_worktree_dirty', False)}`"
        ),
        (
            "- Inventory tool version: "
            f"`{TOOL_VERSION}`"
        ),
        (
            "- Generated: "
            f"`{metadata.get('generated_at_utc', '')}`"
        ),
        "",
        "## Scan summary",
        "",
        f"- Scripts scanned: **{len(rows)}**",
        (
            "- Scripts with detected randomization: "
            f"**{random_count}**"
        ),
        (
            "- Scripts with detected version text: "
            f"**{version_count}**"
        ),
        (
            "- Numeric script IDs represented by "
            "multiple files: "
            f"**{duplicate_number_count}**"
        ),
        "",
        "### Automated axis suggestions",
        "",
    ]

    if axis_counts:
        for axis, count in sorted(
            axis_counts.items()
        ):
            lines.append(
                f"- `{axis}`: {count}"
            )
    else:
        lines.append(
            "- No axis suggestions detected."
        )

    lines.extend(
        [
            "",
            "## Script inventory",
            "",
            (
                "| # | Script | Version | SHA-256 | "
                "Git commit | Random | Automated axis suggestion |"
            ),
            (
                "|---:|---|---|---|---|:---:|---|"
            ),
        ]
    )

    for row in rows:
        columns = [
            markdown_escape(
                row["script_number"]
            ),
            (
                "`"
                + markdown_escape(
                    row["filename"]
                )
                + "`"
            ),
            markdown_escape(
                row["declared_version"]
            ),
            (
                "`"
                + short_sha(row["sha256"])
                + "`"
            ),
            (
                "`"
                + short_sha(
                    row["last_git_commit"]
                )
                + "`"
            ),
            (
                "yes"
                if row[
                    "randomization_detected"
                ] == "yes"
                else ""
            ),
            markdown_escape(
                row[
                    "auto_axis_suggestion"
                ]
            ),
        ]

        lines.append(
            "| "
            + " | ".join(columns)
            + " |"
        )

    lines.extend(
        [
            "",
            "## Review workflow",
            "",
            (
                "`legacy_analysis_inventory_auto.csv` "
                "is regenerated on every scan."
            ),
            "",
            (
                "`legacy_analysis_inventory_review.csv` "
                "preserves manually reviewed fields across "
                "rescans by matching `relative_path`."
            ),
            "",
            (
                "Manual review must determine which scripts "
                "are authoritative, reusable, superseded, "
                "or study-specific."
            ),
            "",
        ]
    )

    AUTO_MD.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def validate_paths() -> None:
    if not SOURCE_REPO_ROOT.exists():
        raise FileNotFoundError(
            "Source repository not found:\n"
            f"{SOURCE_REPO_ROOT}\n\n"
            "Expected the old and new repositories "
            "to be sibling directories."
        )

    if not SOURCE_SCRIPTS_DIR.exists():
        raise FileNotFoundError(
            "Source scripts directory not found:\n"
            f"{SOURCE_SCRIPTS_DIR}"
        )


def main() -> None:
    validate_paths()

    DOCS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    scripts = discover_scripts()

    if not scripts:
        raise RuntimeError(
            "No supported scripts found under:\n"
            f"{SOURCE_SCRIPTS_DIR}"
        )

    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- legacy analysis inventory"
    )
    print("=" * 80)

    print(
        f"Tool version:      {TOOL_VERSION}"
    )
    print(
        f"Source repository: {SOURCE_REPO_ROOT}"
    )
    print(
        f"New repository:    {NEW_REPO_ROOT}"
    )
    print(
        f"Scripts found:     {len(scripts)}"
    )
    print()

    rows: list[dict[str, str]] = []

    for index, path in enumerate(
        scripts,
        start=1,
    ):
        row = scan_script(path)

        rows.append(row)

        print(
            f"[{index:02d}/{len(scripts):02d}] "
            f"{row['filename']} "
            f"sha={short_sha(row['sha256'])} "
            f"version="
            f"{row['declared_version'] or '-'}"
        )

    annotate_duplicate_numbers(rows)

    auto_fields = list(
        rows[0].keys()
    )

    metadata = source_repo_metadata()

    write_csv(
        AUTO_CSV,
        rows,
        auto_fields,
    )

    write_review_csv(
        rows,
        auto_fields,
    )

    write_markdown(
        rows,
        metadata,
    )

    META_JSON.write_text(
        json.dumps(
            metadata,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print()
    print("Inventory written:")

    print(
        "  "
        + str(
            AUTO_CSV.relative_to(
                NEW_REPO_ROOT
            )
        )
    )

    print(
        "  "
        + str(
            REVIEW_CSV.relative_to(
                NEW_REPO_ROOT
            )
        )
    )

    print(
        "  "
        + str(
            AUTO_MD.relative_to(
                NEW_REPO_ROOT
            )
        )
    )

    print(
        "  "
        + str(
            META_JSON.relative_to(
                NEW_REPO_ROOT
            )
        )
    )

    print()
    print(
        "Next step: manually review "
        "authoritative versions and migration targets."
    )

    print("=" * 80)


if __name__ == "__main__":
    main()
