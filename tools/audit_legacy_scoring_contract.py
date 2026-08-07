from __future__ import annotations

import ast
import csv
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


TOOL_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = ROOT.parent / "paper4_sarcoma_dog"

DOCS_DIR = ROOT / "docs"

OUTPUT_MD = (
    DOCS_DIR
    / "legacy_scoring_contract.md"
)

OUTPUT_CSV = (
    DOCS_DIR
    / "legacy_scoring_functions.csv"
)

OUTPUT_JSON = (
    DOCS_DIR
    / "legacy_scoring_contract.json"
)


SOURCE_SCRIPTS = [
    "scripts/21_finalize_canine_transfer_programs_FIXED_V2.py",
    "scripts/22_prepare_human_osteosarcoma_cohorts.py",
    "scripts/23_external_human_validation.py",
    "scripts/25_prepare_gse39055_third_human_cohort.py",
    "scripts/26_validate_gse39055_rfs.py",
    "scripts/28_conservative_module_preservation_audit.py",
    "scripts/31_gse39055_assay_quality_diagnostic_v2.py",
    "scripts/42_score_ammons_single_cell_localization.py",
    "scripts/44_recompute_ammons_six_dog_localization.py",
    "scripts/46_gse239948_external_canine_representation_v2.py",
    "scripts/49_gse239948_blind_de_novo_rediscovery_FIXED_V2.py",
]


FUNCTION_NAME_KEYWORDS = [
    "score",
    "module",
    "program",
    "standard",
    "zscore",
    "z_score",
    "loading",
    "weight",
    "ortholog",
    "mapping",
    "coverage",
]


SIGNAL_PATTERNS = {
    "mean": re.compile(
        r"\.(?:mean|nanmean)\s*\(|"
        r"\bnp\.nanmean\s*\(|"
        r"\bnp\.mean\s*\(",
        re.IGNORECASE,
    ),
    "std": re.compile(
        r"\.(?:std|nanstd)\s*\(|"
        r"\bnp\.nanstd\s*\(|"
        r"\bnp\.std\s*\(",
        re.IGNORECASE,
    ),
    "zscore": re.compile(
        r"\bzscore\b|"
        r"\bz_score\b|"
        r"\bstandardiz",
        re.IGNORECASE,
    ),
    "weights": re.compile(
        r"\bweight(?:s|ed)?\b",
        re.IGNORECASE,
    ),
    "loadings": re.compile(
        r"\bloading(?:s)?\b",
        re.IGNORECASE,
    ),
    "direction": re.compile(
        r"\bdirection\b|"
        r"\bsign(?:ed)?\b",
        re.IGNORECASE,
    ),
    "matrix_product": re.compile(
        r"\.dot\s*\(|"
        r"(?<![@\w])@(?![\w])",
    ),
    "correlation": re.compile(
        r"\bcorr\b|"
        r"\bpearson\b|"
        r"\bspearman\b",
        re.IGNORECASE,
    ),
    "coverage": re.compile(
        r"\bcoverage\b|"
        r"\boverlap\b",
        re.IGNORECASE,
    ),
    "rank": re.compile(
        r"\brank\b",
        re.IGNORECASE,
    ),
}


CSV_LITERAL_PATTERN = re.compile(
    r"""["']([^"']+\.csv)["']""",
    re.IGNORECASE,
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


def read_text(
    path: Path,
) -> str:
    return path.read_text(
        encoding="utf-8",
        errors="replace",
    )


def relevant_function_name(
    name: str,
) -> bool:
    lowered = name.lower()

    return any(
        keyword in lowered
        for keyword in FUNCTION_NAME_KEYWORDS
    )


def detect_signals(
    text: str,
) -> list[str]:
    signals = []

    for label, pattern in (
        SIGNAL_PATTERNS.items()
    ):
        if pattern.search(text):
            signals.append(label)

    return signals


def extract_functions(
    path: Path,
) -> list[dict[str, object]]:
    text = read_text(path)

    try:
        tree = ast.parse(text)
    except SyntaxError:
        return []

    functions = []

    for node in ast.walk(tree):
        if not isinstance(
            node,
            (
                ast.FunctionDef,
                ast.AsyncFunctionDef,
            ),
        ):
            continue

        source = ast.get_source_segment(
            text,
            node,
        )

        if source is None:
            continue

        signals = detect_signals(
            source
        )

        if not (
            relevant_function_name(
                node.name
            )
            or signals
        ):
            continue

        functions.append(
            {
                "script": (
                    path.relative_to(
                        SOURCE_ROOT
                    ).as_posix()
                ),
                "function": node.name,
                "start_line": (
                    node.lineno
                ),
                "end_line": (
                    getattr(
                        node,
                        "end_lineno",
                        node.lineno,
                    )
                ),
                "signals": signals,
                "source": source,
            }
        )

    return sorted(
        functions,
        key=lambda item: (
            int(item["start_line"])
        ),
    )


def extract_relevant_lines(
    path: Path,
) -> list[dict[str, object]]:
    text = read_text(path)

    rows = []

    for line_number, line in enumerate(
        text.splitlines(),
        start=1,
    ):
        stripped = line.strip()

        if not stripped:
            continue

        signals = detect_signals(
            stripped
        )

        if not signals:
            continue

        rows.append(
            {
                "line": line_number,
                "signals": signals,
                "text": stripped,
            }
        )

    return rows


def extract_csv_literals(
    path: Path,
) -> list[str]:
    text = read_text(path)

    return sorted(
        set(
            match.group(1)
            for match in (
                CSV_LITERAL_PATTERN.finditer(
                    text
                )
            )
        )
    )


def find_result_file(
    literal: str,
) -> Path | None:
    basename = Path(
        literal.replace("\\", "/")
    ).name

    candidates = list(
        (
            SOURCE_ROOT
            / "results"
            / "tables"
        ).glob(
            basename
        )
    )

    if len(candidates) == 1:
        return candidates[0]

    return None


def csv_schema(
    path: Path,
) -> dict[str, object]:
    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        reader = csv.reader(
            handle
        )

        try:
            header = next(
                reader
            )
        except StopIteration:
            header = []

        row_count = sum(
            1
            for _ in reader
        )

    return {
        "path": (
            path.relative_to(
                SOURCE_ROOT
            ).as_posix()
        ),
        "columns": header,
        "row_count": row_count,
        "sha256": sha256_file(
            path
        ),
    }


def script_record(
    relative_path: str,
) -> dict[str, object]:
    path = (
        SOURCE_ROOT
        / relative_path
    )

    if not path.exists():
        raise FileNotFoundError(
            f"Missing source script: "
            f"{path}"
        )

    csv_files = []

    for literal in extract_csv_literals(
        path
    ):
        resolved = find_result_file(
            literal
        )

        if resolved is None:
            continue

        csv_files.append(
            csv_schema(
                resolved
            )
        )

    return {
        "script": relative_path,
        "sha256": sha256_file(
            path
        ),
        "functions": (
            extract_functions(
                path
            )
        ),
        "relevant_lines": (
            extract_relevant_lines(
                path
            )
        ),
        "referenced_csv_schemas": (
            csv_files
        ),
    }


def write_function_csv(
    records: list[
        dict[str, object]
    ],
) -> None:
    rows = []

    for record in records:
        for function in record[
            "functions"
        ]:
            rows.append(
                {
                    "script": function[
                        "script"
                    ],
                    "function": function[
                        "function"
                    ],
                    "start_line": function[
                        "start_line"
                    ],
                    "end_line": function[
                        "end_line"
                    ],
                    "signals": ";".join(
                        function[
                            "signals"
                        ]
                    ),
                }
            )

    fieldnames = [
        "script",
        "function",
        "start_line",
        "end_line",
        "signals",
    ]

    with OUTPUT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(
            rows
        )


def fence_code(
    source: str,
) -> list[str]:
    return [
        "```python",
        source.rstrip(),
        "```",
    ]


def write_markdown(
    records: list[
        dict[str, object]
    ],
) -> None:
    lines = [
        "# Legacy Scoring Contract Audit",
        "",
        (
            "> This report extracts scoring-, "
            "mapping-, weighting-, coverage-, "
            "and standardization-related logic "
            "from authoritative Paper 4 scripts."
        ),
        "",
        (
            "It is descriptive only. "
            "No reusable scoring contract is "
            "defined until these legacy semantics "
            "have been reviewed."
        ),
        "",
        f"- Tool version: `{TOOL_VERSION}`",
        (
            "- Generated: "
            f"`{datetime.now(timezone.utc).isoformat()}`"
        ),
        "",
    ]

    for record in records:
        lines.extend(
            [
                "---",
                "",
                (
                    "## `"
                    + str(
                        record["script"]
                    )
                    + "`"
                ),
                "",
                (
                    "- SHA-256: `"
                    + str(
                        record["sha256"]
                    )
                    + "`"
                ),
                "",
                "### Relevant functions",
                "",
            ]
        )

        functions = record[
            "functions"
        ]

        if not functions:
            lines.append(
                "No candidate functions detected."
            )

        for function in functions:
            signals = ", ".join(
                function[
                    "signals"
                ]
            )

            lines.extend(
                [
                    (
                        "#### `"
                        + str(
                            function[
                                "function"
                            ]
                        )
                        + "`"
                    ),
                    "",
                    (
                        "- Lines: "
                        f"{function['start_line']}"
                        "-"
                        f"{function['end_line']}"
                    ),
                    (
                        "- Signals: "
                        f"`{signals}`"
                    ),
                    "",
                    *fence_code(
                        str(
                            function[
                                "source"
                            ]
                        )
                    ),
                    "",
                ]
            )

        lines.extend(
            [
                "### Relevant standalone lines",
                "",
                "| Line | Signals | Code |",
                "|---:|---|---|",
            ]
        )

        for item in record[
            "relevant_lines"
        ]:
            code = (
                str(
                    item["text"]
                )
                .replace(
                    "|",
                    "\\|",
                )
                .replace(
                    "`",
                    "\\`",
                )
            )

            signals = ", ".join(
                item[
                    "signals"
                ]
            )

            lines.append(
                "| "
                f"{item['line']} | "
                f"{signals} | "
                f"`{code}` |"
            )

        lines.extend(
            [
                "",
                "### Referenced result-table schemas",
                "",
            ]
        )

        schemas = record[
            "referenced_csv_schemas"
        ]

        if not schemas:
            lines.append(
                "No resolvable result-table CSVs detected."
            )

        for schema in schemas:
            lines.extend(
                [
                    (
                        "#### `"
                        + str(
                            schema[
                                "path"
                            ]
                        )
                        + "`"
                    ),
                    "",
                    (
                        f"- Rows: "
                        f"{schema['row_count']}"
                    ),
                    (
                        "- SHA-256: `"
                        + str(
                            schema[
                                "sha256"
                            ]
                        )
                        + "`"
                    ),
                    (
                        "- Columns: "
                        + ", ".join(
                            "`"
                            + str(column)
                            + "`"
                            for column in (
                                schema[
                                    "columns"
                                ]
                            )
                        )
                    ),
                    "",
                ]
            )

    OUTPUT_MD.write_text(
        "\n".join(
            lines
        ),
        encoding="utf-8",
    )


def main() -> None:
    DOCS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    records = []

    print(
        "=" * 80
    )

    print(
        "Molecular Transport Audit "
        "- legacy scoring contract audit"
    )

    print(
        "=" * 80
    )

    for index, relative_path in enumerate(
        SOURCE_SCRIPTS,
        start=1,
    ):
        print(
            f"[{index:02d}/"
            f"{len(SOURCE_SCRIPTS):02d}] "
            f"{relative_path}"
        )

        records.append(
            script_record(
                relative_path
            )
        )

    write_function_csv(
        records
    )

    write_markdown(
        records
    )

    payload = {
        "tool_version": (
            TOOL_VERSION
        ),
        "generated_at_utc": (
            datetime.now(
                timezone.utc
            ).isoformat()
        ),
        "source_repository": str(
            SOURCE_ROOT
        ),
        "records": records,
    }

    OUTPUT_JSON.write_text(
        json.dumps(
            payload,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print()

    print(
        f"Written: "
        f"{OUTPUT_MD.relative_to(ROOT)}"
    )

    print(
        f"Written: "
        f"{OUTPUT_CSV.relative_to(ROOT)}"
    )

    print(
        f"Written: "
        f"{OUTPUT_JSON.relative_to(ROOT)}"
    )

    print(
        "=" * 80
    )


if __name__ == "__main__":
    main()
