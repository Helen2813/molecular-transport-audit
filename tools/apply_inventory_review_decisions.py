from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

CONFIG_PATH = (
    ROOT
    / "configs"
    / "legacy_authoritative_versions.json"
)

REVIEW_PATH = (
    ROOT
    / "docs"
    / "legacy_analysis_inventory_review.csv"
)


def load_config() -> dict:
    return json.loads(
        CONFIG_PATH.read_text(
            encoding="utf-8"
        )
    )


def load_inventory() -> tuple[list[dict[str, str]], list[str]]:
    with REVIEW_PATH.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle)

        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    return rows, fieldnames


def append_note(
    existing: str,
    addition: str,
) -> str:
    existing = existing.strip()
    addition = addition.strip()

    if not existing:
        return addition

    if addition in existing:
        return existing

    return f"{existing} | {addition}"


def apply_duplicate_decisions(
    rows: list[dict[str, str]],
    config: dict,
) -> None:
    by_filename = {
        row["filename"]: row
        for row in rows
    }

    for decision in config["decisions"]:
        authoritative = decision[
            "authoritative_file"
        ]

        if authoritative not in by_filename:
            raise RuntimeError(
                f"Authoritative file missing from inventory: "
                f"{authoritative}"
            )

        row = by_filename[authoritative]

        row["authoritative_version"] = "yes"
        row["review_status"] = "authoritative"

        row["notes"] = append_note(
            row.get("notes", ""),
            decision["reason"],
        )

        expected_version = decision.get(
            "expected_version",
            "",
        )

        observed_version = row.get(
            "declared_version",
            "",
        )

        if (
            expected_version
            and observed_version != expected_version
        ):
            raise RuntimeError(
                "Version mismatch for "
                f"{authoritative}: "
                f"expected={expected_version}, "
                f"observed={observed_version}"
            )

        for superseded in decision.get(
            "superseded_files",
            [],
        ):
            if superseded not in by_filename:
                raise RuntimeError(
                    f"Superseded file missing: "
                    f"{superseded}"
                )

            superseded_row = by_filename[
                superseded
            ]

            superseded_row[
                "authoritative_version"
            ] = "no"

            superseded_row[
                "review_status"
            ] = "superseded"

            superseded_row[
                "notes"
            ] = append_note(
                superseded_row.get(
                    "notes",
                    "",
                ),
                (
                    "Superseded by "
                    f"{authoritative}."
                ),
            )


def apply_additional_authoritative(
    rows: list[dict[str, str]],
    config: dict,
) -> None:
    by_filename = {
        row["filename"]: row
        for row in rows
    }

    for decision in config.get(
        "additional_authoritative_scripts",
        [],
    ):
        filename = decision["file"]

        if filename not in by_filename:
            raise RuntimeError(
                f"File missing from inventory: "
                f"{filename}"
            )

        row = by_filename[filename]

        expected = decision.get(
            "expected_version",
            "",
        )

        observed = row.get(
            "declared_version",
            "",
        )

        if expected and observed != expected:
            raise RuntimeError(
                f"Version mismatch for {filename}: "
                f"expected={expected}, "
                f"observed={observed}"
            )

        row[
            "authoritative_version"
        ] = "yes"

        row[
            "review_status"
        ] = "authoritative"

        row["notes"] = append_note(
            row.get("notes", ""),
            decision["reason"],
        )


def write_inventory(
    rows: list[dict[str, str]],
    fieldnames: list[str],
) -> None:
    with REVIEW_PATH.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    config = load_config()

    rows, fieldnames = load_inventory()

    apply_duplicate_decisions(
        rows,
        config,
    )

    apply_additional_authoritative(
        rows,
        config,
    )

    write_inventory(
        rows,
        fieldnames,
    )

    authoritative = [
        row
        for row in rows
        if row.get(
            "review_status"
        ) == "authoritative"
    ]

    superseded = [
        row
        for row in rows
        if row.get(
            "review_status"
        ) == "superseded"
    ]

    print("=" * 80)
    print(
        "Legacy inventory review decisions applied"
    )
    print("=" * 80)

    print(
        f"Authoritative files: "
        f"{len(authoritative)}"
    )

    for row in authoritative:
        print(
            f"  [KEEP] {row['filename']}"
        )

    print()

    print(
        f"Superseded files: "
        f"{len(superseded)}"
    )

    for row in superseded:
        print(
            f"  [DROP] {row['filename']}"
        )

    print("=" * 80)


if __name__ == "__main__":
    main()
