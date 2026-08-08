"""Extend the locked preservation fixture with reliability artifacts."""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


TOOL_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]

SOURCE_ROOT = (
    ROOT.parent
    / "paper4_sarcoma_dog"
)

SOURCE_LOO = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE239948_external_module_gene_leave_one_out.csv"
)

SOURCE_RELIABILITY = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE239948_external_module_score_reliability.csv"
)

SOURCE_MANIFEST = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE239948_external_representation_manifest.json"
)

FIXTURE_DIR = (
    ROOT
    / "reference_results"
    / "osteosarcoma_locked"
    / "preservation_fixture"
)

BASE_FIXTURE_MANIFEST = (
    FIXTURE_DIR
    / "fixture_manifest.json"
)

EXISTING_RELIABILITY = (
    FIXTURE_DIR
    / "expected_split_half_reliability.csv"
)

OUTPUT_LOO = (
    FIXTURE_DIR
    / "expected_gene_leave_one_out.csv"
)

OUTPUT_MANIFEST = (
    FIXTURE_DIR
    / "reliability_fixture_manifest.json"
)

EXPECTED_SCRIPT_VERSION = (
    "46-gse239948-external-canine-representation-v2"
)

PRIMARY_MODULES = [
    "M34",
    "M11",
    "M24",
    "M40",
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


def read_json(
    path: Path,
) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(
            f"Required JSON missing:\n{path}"
        )

    payload = json.loads(
        path.read_text(
            encoding="utf-8"
        )
    )

    if not isinstance(
        payload,
        dict,
    ):
        raise ValueError(
            f"Expected JSON object: {path}"
        )

    return payload


def manifest_hash_by_basename(
    manifest: dict[str, Any],
    basename: str,
) -> str:
    outputs = manifest.get(
        "outputs",
        {},
    )

    if not isinstance(
        outputs,
        dict,
    ):
        raise ValueError(
            "Legacy manifest has no outputs dictionary."
        )

    matches: list[str] = []

    for raw_path, metadata in outputs.items():
        if (
            Path(str(raw_path)).name
            != basename
        ):
            continue

        if not isinstance(
            metadata,
            dict,
        ):
            continue

        value = metadata.get(
            "sha256"
        )

        if value:
            matches.append(
                str(value)
            )

    matches = sorted(
        set(matches)
    )

    if len(matches) != 1:
        raise RuntimeError(
            "Could not resolve exactly one "
            f"manifest hash for {basename}: "
            f"{matches}"
        )

    return matches[0]


def verify_legacy_output(
    manifest: dict[str, Any],
    path: Path,
) -> str:
    expected = (
        manifest_hash_by_basename(
            manifest,
            path.name,
        )
    )

    observed = sha256_file(
        path
    )

    if expected != observed:
        raise RuntimeError(
            "Legacy output hash mismatch:\n"
            f"file={path}\n"
            f"expected={expected}\n"
            f"observed={observed}"
        )

    return observed


def file_record(
    path: Path,
) -> dict[str, Any]:
    table = pd.read_csv(
        path,
        low_memory=False,
    )

    return {
        "sha256":
            sha256_file(path),
        "size_bytes":
            path.stat().st_size,
        "rows":
            int(table.shape[0]),
        "columns":
            int(table.shape[1]),
    }


def main() -> None:
    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- build preservation reliability fixture"
    )
    print("=" * 80)

    for path in [
        SOURCE_LOO,
        SOURCE_RELIABILITY,
        SOURCE_MANIFEST,
        BASE_FIXTURE_MANIFEST,
        EXISTING_RELIABILITY,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Required file missing:\n{path}"
            )

    legacy_manifest = (
        read_json(
            SOURCE_MANIFEST
        )
    )

    script_version = str(
        legacy_manifest.get(
            "script_version",
            "",
        )
    )

    if (
        script_version
        != EXPECTED_SCRIPT_VERSION
    ):
        raise RuntimeError(
            "Unexpected legacy script version:\n"
            f"expected={EXPECTED_SCRIPT_VERSION}\n"
            f"observed={script_version}"
        )

    if bool(
        legacy_manifest.get(
            "outcome_loaded",
            True,
        )
    ):
        raise RuntimeError(
            "Legacy manifest does not confirm "
            "outcome_loaded=false."
        )

    if int(
        legacy_manifest.get(
            "split_half_repeats",
            -1,
        )
    ) != 2000:
        raise RuntimeError(
            "Legacy manifest does not confirm "
            "2000 split-half repeats."
        )

    print("")
    print(
        "Verifying authoritative reliability "
        "outputs against script 46 manifest..."
    )

    reliability_sha = (
        verify_legacy_output(
            legacy_manifest,
            SOURCE_RELIABILITY,
        )
    )

    loo_sha = verify_legacy_output(
        legacy_manifest,
        SOURCE_LOO,
    )

    print(
        "Legacy reliability hashes: PASS"
    )

    source_reliability = pd.read_csv(
        SOURCE_RELIABILITY
    )

    fixture_reliability = pd.read_csv(
        EXISTING_RELIABILITY
    )

    if (
        list(source_reliability.columns)
        != list(
            fixture_reliability.columns
        )
    ):
        raise RuntimeError(
            "Existing fixture reliability "
            "schema differs from legacy."
        )

    source_ordered = (
        source_reliability
        .sort_values(
            "module_label"
        )
        .reset_index(
            drop=True
        )
    )

    fixture_ordered = (
        fixture_reliability
        .sort_values(
            "module_label"
        )
        .reset_index(
            drop=True
        )
    )

    if not source_ordered.equals(
        fixture_ordered
    ):
        raise RuntimeError(
            "Existing split-half fixture does "
            "not exactly equal the authoritative "
            "legacy reliability table."
        )

    loo = pd.read_csv(
        SOURCE_LOO
    )

    required_loo = {
        "module_label",
        "left_out_gene",
        "n_genes_remaining",
        "correlation_with_full_score",
    }

    missing = sorted(
        required_loo.difference(
            loo.columns
        )
    )

    if missing:
        raise ValueError(
            "Legacy LOO table is missing: "
            + ", ".join(missing)
        )

    modules = set(
        loo[
            "module_label"
        ].astype(str)
    )

    unexpected_modules = sorted(
        modules.difference(
            PRIMARY_MODULES
        )
    )

    if unexpected_modules:
        raise RuntimeError(
            "Unexpected modules in LOO table: "
            + ", ".join(
                unexpected_modules
            )
        )

    shutil.copy2(
        SOURCE_LOO,
        OUTPUT_LOO,
    )

    if (
        sha256_file(
            OUTPUT_LOO
        )
        != loo_sha
    ):
        raise RuntimeError(
            "Copied LOO fixture hash mismatch."
        )

    reliability = (
        fixture_reliability
        .set_index(
            "module_label"
        )
        .reindex(
            PRIMARY_MODULES
        )
        .reset_index()
    )

    print("")
    print(
        "Locked reliability summary:"
    )
    print("")

    print(
        reliability[
            [
                "module_label",
                "n_common_genes",
                "split_half_median",
                "split_half_q05",
                "split_half_q95",
                "split_half_valid_repeats",
                "minimum_gene_loo_correlation",
                "median_gene_loo_correlation",
            ]
        ].to_string(
            index=False
        )
    )

    loo_counts = (
        loo.groupby(
            "module_label"
        )
        .size()
        .reindex(
            PRIMARY_MODULES
        )
    )

    print("")
    print(
        "Gene leave-one-out rows:"
    )
    print("")

    print(
        loo_counts.to_string()
    )

    manifest = {
        "fixture_builder_version":
            TOOL_VERSION,
        "created_at_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),
        "source_script_version":
            EXPECTED_SCRIPT_VERSION,
        "source_legacy_manifest_sha256":
            sha256_file(
                SOURCE_MANIFEST
            ),
        "base_preservation_fixture_manifest_sha256":
            sha256_file(
                BASE_FIXTURE_MANIFEST
            ),
        "outcome_loaded":
            False,
        "contract": {
            "modules":
                PRIMARY_MODULES,
            "split_half_repeats":
                2000,
            "base_seed":
                42,
            "module_seed_stride":
                100,
            "split_half_seed_offset":
                2,
            "split_half_seeds": {
                "M34": 44,
                "M11": 144,
                "M24": 244,
                "M40": 344,
            },
            "split_rule":
                (
                    "rng_permutation_of_gene_order_"
                    "then_floor_half_nonoverlapping"
                ),
            "half_score":
                "frozen_signed_mean_unstandardized",
            "half_agreement":
                "pearson_correlation",
            "summary":
                "median_q05_q95",
            "gene_loo":
                (
                    "pearson_correlation_of_each_"
                    "leave_one_gene_out_score_with_"
                    "full_frozen_signed_score"
                ),
            "minimum_module_genes":
                3,
        },
        "authoritative_source_outputs": {
            SOURCE_RELIABILITY.name: {
                "sha256":
                    reliability_sha,
            },
            SOURCE_LOO.name: {
                "sha256":
                    loo_sha,
            },
        },
        "fixture_files": {
            EXISTING_RELIABILITY.name:
                file_record(
                    EXISTING_RELIABILITY
                ),
            OUTPUT_LOO.name:
                file_record(
                    OUTPUT_LOO
                ),
        },
    }

    OUTPUT_MANIFEST.write_text(
        json.dumps(
            manifest,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("")
    print("=" * 80)
    print(
        "Reliability fixture summary"
    )
    print("=" * 80)
    print(
        "Split-half repeats: 2000"
    )
    print(
        "Seeds: "
        "M34=44, M11=144, "
        "M24=244, M40=344"
    )
    print(
        f"Total LOO rows: "
        f"{loo.shape[0]}"
    )
    print(
        "Manifest: "
        + str(
            OUTPUT_MANIFEST.relative_to(
                ROOT
            )
        )
    )
    print("=" * 80)


if __name__ == "__main__":
    main()
