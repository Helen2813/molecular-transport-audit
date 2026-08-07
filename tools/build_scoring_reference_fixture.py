from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


TOOL_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = ROOT.parent / "paper4_sarcoma_dog"

RESULTS_DIR = SOURCE_ROOT / "results" / "tables"
HUMAN_DIR = (
    SOURCE_ROOT
    / "data"
    / "processed"
    / "human_validation"
)

STRICT_WEIGHTS_FILE = (
    RESULTS_DIR
    / "GSE238110_frozen_transfer_gene_weights_strict.csv"
)

SCORING_SPEC_FILE = (
    RESULTS_DIR
    / "GSE238110_frozen_transfer_scoring_specification.csv"
)

PREPARATION_MANIFEST_FILE = (
    RESULTS_DIR
    / "human_validation_cohort_preparation_manifest.json"
)

TARGET_EXPRESSION_FILE = (
    HUMAN_DIR
    / "TARGET_OS_expression_log2_gene_symbol.csv"
)

TARGET_SCORES_FILE = (
    HUMAN_DIR
    / "TARGET_OS_frozen_transfer_scores.csv"
)

TARGET_COVERAGE_FILE = (
    RESULTS_DIR
    / "TARGET_OS_frozen_transfer_score_coverage.csv"
)

OUTPUT_DIR = (
    ROOT
    / "reference_results"
    / "osteosarcoma_locked"
    / "scoring_fixture"
)

PRIMARY_MODULES = [
    "M34",
    "M11",
    "M24",
    "M40",
]

EXPECTED_SCORE_SUFFIXES = [
    "__strict__signed_mean_z",
    "__strict__canine_pca_weighted_z",
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


def run_git(args: list[str]) -> str:
    completed = subprocess.run(
        ["git", *args],
        cwd=SOURCE_ROOT,
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )

    return completed.stdout.strip()


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"Required source file not found: {path}"
        )


def verify_prepared_file(
    path: Path,
    preparation_manifest: dict,
) -> dict[str, str]:
    files = preparation_manifest.get(
        "files",
        {},
    )

    observed = sha256_file(path)

    expected = (
        files.get(
            path.name,
            {},
        ).get(
            "sha256",
            ""
        )
    )

    if expected and observed != expected:
        raise RuntimeError(
            "Prepared source hash mismatch:\n"
            f"file={path}\n"
            f"expected={expected}\n"
            f"observed={observed}"
        )

    return {
        "path": str(
            path.relative_to(
                SOURCE_ROOT
            )
        ),
        "sha256": observed,
        "manifest_expected_sha256": expected,
        "manifest_hash_verified": bool(
            expected
            and expected == observed
        ),
    }


def write_csv(
    table: pd.DataFrame,
    path: Path,
    include_index: bool,
) -> None:
    table.to_csv(
        path,
        index=include_index,
        index_label=(
            table.index.name
            or "sample_id"
            if include_index
            else None
        ),
    )


def main() -> None:
    required = [
        STRICT_WEIGHTS_FILE,
        SCORING_SPEC_FILE,
        PREPARATION_MANIFEST_FILE,
        TARGET_EXPRESSION_FILE,
        TARGET_SCORES_FILE,
        TARGET_COVERAGE_FILE,
    ]

    for path in required:
        require_file(path)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    preparation_manifest = json.loads(
        PREPARATION_MANIFEST_FILE.read_text(
            encoding="utf-8"
        )
    )

    source_expression_record = (
        verify_prepared_file(
            TARGET_EXPRESSION_FILE,
            preparation_manifest,
        )
    )

    source_scores_record = (
        verify_prepared_file(
            TARGET_SCORES_FILE,
            preparation_manifest,
        )
    )

    strict_weights = pd.read_csv(
        STRICT_WEIGHTS_FILE
    )

    scoring_spec = pd.read_csv(
        SCORING_SPEC_FILE
    )

    expression = pd.read_csv(
        TARGET_EXPRESSION_FILE,
        index_col=0,
    )

    expected_scores = pd.read_csv(
        TARGET_SCORES_FILE,
        index_col=0,
    )

    coverage = pd.read_csv(
        TARGET_COVERAGE_FILE
    )

    primary_weights = strict_weights[
        strict_weights[
            "module_label"
        ].isin(
            PRIMARY_MODULES
        )
    ].copy()

    if primary_weights.empty:
        raise RuntimeError(
            "No primary frozen strict weights found."
        )

    primary_weights[
        "human_gene_symbol"
    ] = (
        primary_weights[
            "human_gene_symbol"
        ]
        .astype(str)
        .str.upper()
    )

    requested_genes = (
        primary_weights[
            "human_gene_symbol"
        ]
        .drop_duplicates()
        .tolist()
    )

    expression = expression.copy()

    expression.columns = (
        expression.columns
        .astype(str)
        .str.upper()
    )

    expression = expression.loc[
        :,
        ~expression.columns.duplicated(
            keep="first"
        ),
    ]

    available_genes = [
        gene
        for gene in requested_genes
        if gene in expression.columns
    ]

    missing_genes = [
        gene
        for gene in requested_genes
        if gene not in expression.columns
    ]

    if not available_genes:
        raise RuntimeError(
            "No primary frozen genes are present "
            "in the TARGET expression matrix."
        )

    fixture_expression = (
        expression[
            available_genes
        ].copy()
    )

    score_columns = []

    for module in PRIMARY_MODULES:
        for suffix in EXPECTED_SCORE_SUFFIXES:
            column = (
                module
                + suffix
            )

            if column not in expected_scores.columns:
                raise RuntimeError(
                    "Expected legacy score column "
                    f"is missing: {column}"
                )

            score_columns.append(
                column
            )

    common_samples = (
        fixture_expression.index
        .intersection(
            expected_scores.index
        )
    )

    if common_samples.empty:
        raise RuntimeError(
            "Expression and legacy score tables "
            "have no shared TARGET samples."
        )

    fixture_expression = (
        fixture_expression.loc[
            common_samples
        ].copy()
    )

    fixture_scores = (
        expected_scores.loc[
            common_samples,
            score_columns,
        ].copy()
    )

    fixture_spec = scoring_spec[
        scoring_spec[
            "module_label"
        ].isin(
            PRIMARY_MODULES
        )
    ].copy()

    fixture_coverage = coverage[
        coverage[
            "module_label"
        ].isin(
            PRIMARY_MODULES
        )
    ].copy()

    if (
        "mapping" in
        fixture_coverage.columns
    ):
        fixture_coverage = (
            fixture_coverage[
                fixture_coverage[
                    "mapping"
                ].astype(str)
                .str.lower()
                .eq("strict")
            ].copy()
        )

    expression_output = (
        OUTPUT_DIR
        / "TARGET_OS_primary_strict_expression.csv"
    )

    weights_output = (
        OUTPUT_DIR
        / "primary_strict_gene_weights.csv"
    )

    scores_output = (
        OUTPUT_DIR
        / "TARGET_OS_expected_primary_scores.csv"
    )

    coverage_output = (
        OUTPUT_DIR
        / "TARGET_OS_expected_strict_coverage.csv"
    )

    spec_output = (
        OUTPUT_DIR
        / "primary_scoring_specification.csv"
    )

    manifest_output = (
        OUTPUT_DIR
        / "fixture_manifest.json"
    )

    write_csv(
        fixture_expression,
        expression_output,
        include_index=True,
    )

    write_csv(
        primary_weights,
        weights_output,
        include_index=False,
    )

    write_csv(
        fixture_scores,
        scores_output,
        include_index=True,
    )

    write_csv(
        fixture_coverage,
        coverage_output,
        include_index=False,
    )

    write_csv(
        fixture_spec,
        spec_output,
        include_index=False,
    )

    source_head = run_git(
        [
            "rev-parse",
            "HEAD",
        ]
    )

    manifest = {
        "fixture_version": "0.1.0",
        "builder_version": TOOL_VERSION,
        "generated_at_utc": (
            datetime.now(
                timezone.utc
            ).isoformat()
        ),
        "source_repository": (
            SOURCE_ROOT.name
        ),
        "source_repository_head": (
            source_head
        ),
        "cohort": "TARGET_OS",
        "primary_modules": (
            PRIMARY_MODULES
        ),
        "scoring_contract": {
            "primary_score": (
                "strict_one_to_one_signed_mean_z"
            ),
            "secondary_score": (
                "strict_one_to_one_canine_pca_weighted_z"
            ),
            "minimum_genes": 3,
            "minimum_fraction": 0.50,
            "gene_standardization": (
                "within_cohort_gene_zscore_ddof1"
            ),
            "score_standardization": (
                "within_cohort_score_zscore_ddof1"
            ),
        },
        "source_files": {
            "expression": (
                source_expression_record
            ),
            "legacy_scores": (
                source_scores_record
            ),
            "strict_weights": {
                "path": str(
                    STRICT_WEIGHTS_FILE.relative_to(
                        SOURCE_ROOT
                    )
                ),
                "sha256": sha256_file(
                    STRICT_WEIGHTS_FILE
                ),
            },
            "scoring_specification": {
                "path": str(
                    SCORING_SPEC_FILE.relative_to(
                        SOURCE_ROOT
                    )
                ),
                "sha256": sha256_file(
                    SCORING_SPEC_FILE
                ),
            },
            "preparation_manifest": {
                "path": str(
                    PREPARATION_MANIFEST_FILE.relative_to(
                        SOURCE_ROOT
                    )
                ),
                "sha256": sha256_file(
                    PREPARATION_MANIFEST_FILE
                ),
            },
        },
        "fixture_summary": {
            "n_samples": int(
                fixture_expression.shape[0]
            ),
            "n_expression_genes": int(
                fixture_expression.shape[1]
            ),
            "n_weight_rows": int(
                primary_weights.shape[0]
            ),
            "n_expected_score_columns": int(
                fixture_scores.shape[1]
            ),
            "n_missing_requested_genes": int(
                len(missing_genes)
            ),
            "missing_requested_genes": (
                missing_genes
            ),
        },
        "fixture_files": {},
    }

    fixture_paths = [
        expression_output,
        weights_output,
        scores_output,
        coverage_output,
        spec_output,
    ]

    for path in fixture_paths:
        manifest[
            "fixture_files"
        ][path.name] = {
            "sha256": sha256_file(
                path
            ),
            "size_bytes": (
                path.stat().st_size
            ),
        }

    manifest_output.write_text(
        json.dumps(
            manifest,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("=" * 80)

    print(
        "Molecular Transport Audit "
        "- TARGET scoring regression fixture"
    )

    print("=" * 80)

    print(
        f"Source HEAD: {source_head}"
    )

    print(
        f"Samples: {fixture_expression.shape[0]}"
    )

    print(
        "Primary strict genes available: "
        f"{fixture_expression.shape[1]}"
    )

    print(
        "Primary weight rows: "
        f"{primary_weights.shape[0]}"
    )

    print(
        "Expected score columns: "
        f"{fixture_scores.shape[1]}"
    )

    print(
        "Missing requested genes: "
        f"{len(missing_genes)}"
    )

    print()

    print(
        f"Written: "
        f"{OUTPUT_DIR.relative_to(ROOT)}"
    )

    print("=" * 80)


if __name__ == "__main__":
    main()
