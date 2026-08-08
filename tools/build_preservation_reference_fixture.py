"""Build the locked GSE239948 preservation regression fixture."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


TOOL_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]

SOURCE_ROOT = (
    ROOT.parent
    / "paper4_sarcoma_dog"
)

SOURCE_SCRIPT = (
    SOURCE_ROOT
    / "scripts"
    / "46_gse239948_external_canine_representation_v2.py"
)

REFERENCE_EXPRESSION_FILE = (
    SOURCE_ROOT
    / "data"
    / "processed"
    / "GSE238110_DOG2_expression_log2cpm_matched_allgenes.csv"
)

EXTERNAL_EXPRESSION_FILE = (
    SOURCE_ROOT
    / "data"
    / "processed"
    / "canine_validation_GSE239948_expression_log2_symbol.csv"
)

STRICT_WEIGHTS_FILE = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE238110_frozen_transfer_gene_weights_strict.csv"
)

GENE_MAP_FILE = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE239948_external_gene_mapping.csv"
)

EXPECTED_COVERAGE_FILE = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE239948_external_frozen_program_coverage.csv"
)

EXPECTED_STRUCTURE_FILE = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE239948_external_module_structure_preservation.csv"
)

EXPECTED_RELIABILITY_FILE = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE239948_external_module_score_reliability.csv"
)

LEGACY_MANIFEST_FILE = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE239948_external_representation_manifest.json"
)

OUTPUT_DIR = (
    ROOT
    / "reference_results"
    / "osteosarcoma_locked"
    / "preservation_fixture"
)

OUTPUT_REFERENCE = (
    OUTPUT_DIR
    / "DOG2_reference_module_expression.csv"
)

OUTPUT_EXTERNAL = (
    OUTPUT_DIR
    / "GSE239948_external_module_expression.csv"
)

OUTPUT_WEIGHTS = (
    OUTPUT_DIR
    / "primary_canine_program_weights.csv"
)

OUTPUT_COVERAGE = (
    OUTPUT_DIR
    / "expected_preservation_coverage.csv"
)

OUTPUT_STRUCTURE = (
    OUTPUT_DIR
    / "expected_direct_preservation.csv"
)

OUTPUT_RELIABILITY = (
    OUTPUT_DIR
    / "expected_split_half_reliability.csv"
)

OUTPUT_MANIFEST = (
    OUTPUT_DIR
    / "fixture_manifest.json"
)

PRIMARY_MODULES = [
    "M34",
    "M11",
    "M24",
    "M40",
]

EXPECTED_SCRIPT_VERSION = (
    "46-gse239948-external-canine-representation-v2"
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


def read_required_csv(
    path: Path,
    *,
    index_col: int | None = None,
) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Required file not found:\n"
            f"{path}"
        )

    return pd.read_csv(
        path,
        index_col=index_col,
        low_memory=False,
    )


def read_required_json(
    path: Path,
) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(
            f"Required JSON not found:\n"
            f"{path}"
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
            f"Expected JSON object in "
            f"{path}"
        )

    return payload


def require_columns(
    table: pd.DataFrame,
    columns: list[str],
    label: str,
) -> None:
    missing = [
        column
        for column in columns
        if column not in table.columns
    ]

    if missing:
        raise ValueError(
            f"{label} is missing columns: "
            + ", ".join(missing)
        )


def run_git(
    args: list[str],
) -> str:
    try:
        completed = subprocess.run(
            [
                "git",
                *args,
            ],
            cwd=SOURCE_ROOT,
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
    except (
        FileNotFoundError,
        subprocess.CalledProcessError,
    ):
        return ""

    return completed.stdout.strip()


def parse_gene_list(
    value: Any,
) -> list[str]:
    if pd.isna(value):
        return []

    result: list[str] = []

    for raw in str(value).split(";"):
        gene = (
            raw.strip()
            .upper()
        )

        if (
            gene
            and gene
            not in result
        ):
            result.append(gene)

    return result


def legacy_manifest_hash(
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
            "Legacy manifest has no "
            "outputs dictionary."
        )

    matches: list[str] = []

    for raw_path, payload in outputs.items():
        if (
            Path(
                str(raw_path)
            ).name
            != basename
        ):
            continue

        if not isinstance(
            payload,
            dict,
        ):
            continue

        value = payload.get(
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
            "legacy manifest hash for "
            f"{basename}. "
            f"Found: {matches}"
        )

    return matches[0]


def verify_legacy_output(
    manifest: dict[str, Any],
    path: Path,
) -> dict[str, Any]:
    expected = (
        legacy_manifest_hash(
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

    return {
        "path":
            str(
                path.relative_to(
                    SOURCE_ROOT
                )
            ),
        "sha256":
            observed,
        "size_bytes":
            path.stat().st_size,
        "legacy_manifest_verified":
            True,
    }


def ordered_common_gene_union(
    coverage: pd.DataFrame,
) -> list[str]:
    require_columns(
        coverage,
        [
            "module_label",
            "common_genes",
        ],
        "Legacy coverage",
    )

    genes: list[str] = []

    indexed = coverage.set_index(
        "module_label"
    )

    for module in PRIMARY_MODULES:
        if module not in indexed.index:
            raise ValueError(
                "Legacy coverage is missing "
                f"module {module}."
            )

        module_genes = parse_gene_list(
            indexed.loc[
                module,
                "common_genes",
            ]
        )

        for gene in module_genes:
            if gene not in genes:
                genes.append(gene)

    if not genes:
        raise RuntimeError(
            "No common frozen genes were "
            "found in the legacy coverage."
        )

    return genes


def prepare_reference_expression(
    reference: pd.DataFrame,
    gene_map: pd.DataFrame,
    genes: list[str],
) -> pd.DataFrame:
    require_columns(
        gene_map,
        [
            "analysis_gene_id",
            "reference_expression_column",
        ],
        "Legacy gene map",
    )

    mapping = gene_map[
        [
            "analysis_gene_id",
            "reference_expression_column",
        ]
    ].copy()

    mapping[
        "analysis_gene_id"
    ] = (
        mapping[
            "analysis_gene_id"
        ]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    mapping[
        "reference_expression_column"
    ] = (
        mapping[
            "reference_expression_column"
        ]
        .astype(str)
        .str.strip()
    )

    mapping = mapping[
        mapping[
            "analysis_gene_id"
        ].isin(genes)
    ].copy()

    mapping = mapping[
        ~mapping[
            "reference_expression_column"
        ].isin(
            [
                "",
                "nan",
                "None",
            ]
        )
    ].copy()

    mapping = (
        mapping.drop_duplicates(
            "analysis_gene_id",
            keep="first",
        )
        .set_index(
            "analysis_gene_id"
        )
    )

    output = pd.DataFrame(
        index=reference.index
    )

    missing: list[str] = []

    reference_columns = set(
        reference.columns.astype(str)
    )

    for gene in genes:
        source_column = ""

        if gene in mapping.index:
            candidate = str(
                mapping.loc[
                    gene,
                    "reference_expression_column",
                ]
            )

            if candidate in reference_columns:
                source_column = candidate

        if (
            not source_column
            and gene
            in reference_columns
        ):
            source_column = gene

        if not source_column:
            missing.append(gene)
            continue

        output[gene] = pd.to_numeric(
            reference[
                source_column
            ],
            errors="coerce",
        )

    if missing:
        raise RuntimeError(
            "Could not reconstruct DOG2 "
            "reference expression for "
            f"{len(missing)} locked common genes:\n"
            + ", ".join(missing)
        )

    return output


def prepare_external_expression(
    external: pd.DataFrame,
    genes: list[str],
) -> pd.DataFrame:
    result = external.copy()

    result.columns = (
        result.columns
        .astype(str)
        .str.strip()
        .str.upper()
    )

    result = result.loc[
        :,
        ~pd.Index(
            result.columns
        ).duplicated(
            keep="first"
        ),
    ].copy()

    missing = [
        gene
        for gene in genes
        if gene
        not in result.columns
    ]

    if missing:
        raise RuntimeError(
            "External processed matrix is "
            "missing locked common genes:\n"
            + ", ".join(missing)
        )

    result = result.loc[
        :,
        genes,
    ].apply(
        pd.to_numeric,
        errors="coerce",
    )

    return result


def audit_fixture_variability(
    reference: pd.DataFrame,
    external: pd.DataFrame,
) -> None:
    reference_variable = (
        reference.var(
            axis=0,
            ddof=1,
        )
        .replace(
            [np.inf, -np.inf],
            np.nan,
        )
    )

    external_variable = (
        external.var(
            axis=0,
            ddof=1,
        )
        .replace(
            [np.inf, -np.inf],
            np.nan,
        )
    )

    bad_reference = (
        reference_variable[
            reference_variable.isna()
            | reference_variable.le(0)
        ]
        .index
        .tolist()
    )

    bad_external = (
        external_variable[
            external_variable.isna()
            | external_variable.le(0)
        ]
        .index
        .tolist()
    )

    if (
        bad_reference
        or bad_external
    ):
        raise RuntimeError(
            "Fixture contains genes that "
            "would be removed by legacy "
            "variance filtering.\n"
            f"Reference: {bad_reference}\n"
            f"External: {bad_external}"
        )


def filter_module_rows(
    table: pd.DataFrame,
) -> pd.DataFrame:
    require_columns(
        table,
        ["module_label"],
        "Module table",
    )

    result = table[
        table[
            "module_label"
        ]
        .astype(str)
        .isin(
            PRIMARY_MODULES
        )
    ].copy()

    order = {
        module: index
        for index, module
        in enumerate(
            PRIMARY_MODULES
        )
    }

    result[
        "_module_order"
    ] = (
        result[
            "module_label"
        ]
        .astype(str)
        .map(order)
    )

    result = (
        result.sort_values(
            "_module_order"
        )
        .drop(
            columns=[
                "_module_order"
            ]
        )
        .reset_index(
            drop=True
        )
    )

    return result


def write_fixture_file_info(
    path: Path,
) -> dict[str, Any]:
    table = pd.read_csv(
        path,
        low_memory=False,
    )

    return {
        "sha256":
            sha256_file(
                path
            ),
        "size_bytes":
            path.stat().st_size,
        "rows":
            int(
                table.shape[0]
            ),
        "columns":
            int(
                table.shape[1]
            ),
    }


def main() -> None:
    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- build preservation reference fixture"
    )
    print("=" * 80)
    print(
        f"Tool version: "
        f"{TOOL_VERSION}"
    )
    print(
        f"Source repository: "
        f"{SOURCE_ROOT}"
    )

    for path in [
        SOURCE_SCRIPT,
        REFERENCE_EXPRESSION_FILE,
        EXTERNAL_EXPRESSION_FILE,
        STRICT_WEIGHTS_FILE,
        GENE_MAP_FILE,
        EXPECTED_COVERAGE_FILE,
        EXPECTED_STRUCTURE_FILE,
        EXPECTED_RELIABILITY_FILE,
        LEGACY_MANIFEST_FILE,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Required source file "
                f"not found:\n{path}"
            )

    legacy_manifest = (
        read_required_json(
            LEGACY_MANIFEST_FILE
        )
    )

    observed_version = str(
        legacy_manifest.get(
            "script_version",
            "",
        )
    )

    if (
        observed_version
        != EXPECTED_SCRIPT_VERSION
    ):
        raise RuntimeError(
            "Unexpected legacy preservation "
            "script version:\n"
            f"expected="
            f"{EXPECTED_SCRIPT_VERSION}\n"
            f"observed="
            f"{observed_version}"
        )

    if bool(
        legacy_manifest.get(
            "outcome_loaded",
            True,
        )
    ):
        raise RuntimeError(
            "Legacy preservation manifest "
            "does not confirm "
            "outcome_loaded=false."
        )

    print("")
    print(
        "Verifying legacy script 46 "
        "outputs against its manifest..."
    )

    verified_outputs = {
        path.name:
            verify_legacy_output(
                legacy_manifest,
                path,
            )
        for path in [
            EXTERNAL_EXPRESSION_FILE,
            EXPECTED_COVERAGE_FILE,
            EXPECTED_STRUCTURE_FILE,
            EXPECTED_RELIABILITY_FILE,
        ]
    }

    print(
        "Legacy output hashes: PASS"
    )

    reference = (
        read_required_csv(
            REFERENCE_EXPRESSION_FILE,
            index_col=0,
        )
    )

    external = (
        read_required_csv(
            EXTERNAL_EXPRESSION_FILE,
            index_col=0,
        )
    )

    weights = (
        read_required_csv(
            STRICT_WEIGHTS_FILE
        )
    )

    gene_map = (
        read_required_csv(
            GENE_MAP_FILE
        )
    )

    coverage = filter_module_rows(
        read_required_csv(
            EXPECTED_COVERAGE_FILE
        )
    )

    structure = filter_module_rows(
        read_required_csv(
            EXPECTED_STRUCTURE_FILE
        )
    )

    reliability = filter_module_rows(
        read_required_csv(
            EXPECTED_RELIABILITY_FILE
        )
    )

    require_columns(
        weights,
        [
            "module_label",
            "canine_gene_symbol",
            "risk_oriented_loading",
        ],
        "Frozen strict weights",
    )

    primary_weights = weights[
        weights[
            "module_label"
        ]
        .astype(str)
        .isin(
            PRIMARY_MODULES
        )
    ].copy()

    primary_weights[
        "canine_gene_symbol"
    ] = (
        primary_weights[
            "canine_gene_symbol"
        ]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    genes = (
        ordered_common_gene_union(
            coverage
        )
    )

    reference_fixture = (
        prepare_reference_expression(
            reference,
            gene_map,
            genes,
        )
    )

    external_fixture = (
        prepare_external_expression(
            external,
            genes,
        )
    )

    reference_fixture.index = (
        reference_fixture.index
        .astype(str)
    )

    external_fixture.index = (
        external_fixture.index
        .astype(str)
    )

    audit_fixture_variability(
        reference_fixture,
        external_fixture,
    )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    reference_fixture.to_csv(
        OUTPUT_REFERENCE,
        index=True,
        index_label=(
            reference_fixture.index.name
            or "sample_id"
        ),
    )

    external_fixture.to_csv(
        OUTPUT_EXTERNAL,
        index=True,
        index_label=(
            external_fixture.index.name
            or "sample_id"
        ),
    )

    primary_weights.to_csv(
        OUTPUT_WEIGHTS,
        index=False,
    )

    coverage.to_csv(
        OUTPUT_COVERAGE,
        index=False,
    )

    structure.to_csv(
        OUTPUT_STRUCTURE,
        index=False,
    )

    reliability.to_csv(
        OUTPUT_RELIABILITY,
        index=False,
    )

    expected_metrics_columns = [
        column
        for column in [
            "module_label",
            "n_common_genes",
            "edge_spearman",
            "edge_permutation_p",
            "loading_spearman",
            "loading_permutation_p",
            "external_pc1_variance_explained",
            (
                "pc1_orientation_correlation_"
                "with_frozen_score"
            ),
            "split_half_median",
            "split_half_q05",
            "split_half_q95",
            "split_half_valid_repeats",
            "estimable",
            "nonestimable_reason",
        ]
        if column
        in structure.columns
    ]

    print("")
    print(
        "Locked deterministic/stochastic "
        "preservation targets:"
    )
    print("")

    print(
        structure[
            expected_metrics_columns
        ].to_string(
            index=False
        )
    )

    module_counts = {}

    coverage_index = (
        coverage.set_index(
            "module_label"
        )
    )

    for module in PRIMARY_MODULES:
        module_counts[module] = {
            "n_frozen_genes":
                int(
                    coverage_index.loc[
                        module,
                        "n_frozen_genes",
                    ]
                ),
            "n_common_genes":
                int(
                    coverage_index.loc[
                        module,
                        "n_common_genes",
                    ]
                ),
            "coverage_fraction":
                float(
                    coverage_index.loc[
                        module,
                        "coverage_fraction",
                    ]
                ),
        }

    fixture_files = {
        OUTPUT_REFERENCE.name:
            write_fixture_file_info(
                OUTPUT_REFERENCE
            ),
        OUTPUT_EXTERNAL.name:
            write_fixture_file_info(
                OUTPUT_EXTERNAL
            ),
        OUTPUT_WEIGHTS.name:
            write_fixture_file_info(
                OUTPUT_WEIGHTS
            ),
        OUTPUT_COVERAGE.name:
            write_fixture_file_info(
                OUTPUT_COVERAGE
            ),
        OUTPUT_STRUCTURE.name:
            write_fixture_file_info(
                OUTPUT_STRUCTURE
            ),
        OUTPUT_RELIABILITY.name:
            write_fixture_file_info(
                OUTPUT_RELIABILITY
            ),
    }

    source_head = run_git(
        [
            "rev-parse",
            "HEAD",
        ]
    )

    manifest = {
        "fixture_builder_version":
            TOOL_VERSION,
        "created_at_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),
        "study":
            "canine-human-osteosarcoma",
        "reference_cohort":
            "GSE238110_DOG2",
        "external_cohort":
            "GSE239948",
        "source_repository":
            SOURCE_ROOT.name,
        "source_repository_head":
            source_head,
        "source_script": {
            "path":
                str(
                    SOURCE_SCRIPT.relative_to(
                        SOURCE_ROOT
                    )
                ),
            "script_version":
                EXPECTED_SCRIPT_VERSION,
            "sha256":
                sha256_file(
                    SOURCE_SCRIPT
                ),
        },
        "legacy_manifest": {
            "path":
                str(
                    LEGACY_MANIFEST_FILE
                    .relative_to(
                        SOURCE_ROOT
                    )
                ),
            "sha256":
                sha256_file(
                    LEGACY_MANIFEST_FILE
                ),
            "outcome_loaded":
                False,
            "gene_label_permutations":
                legacy_manifest.get(
                    "gene_label_permutations"
                ),
            "random_panels":
                legacy_manifest.get(
                    "random_panels"
                ),
            "split_half_repeats":
                legacy_manifest.get(
                    "split_half_repeats"
                ),
        },
        "verified_legacy_outputs":
            verified_outputs,
        "source_inputs": {
            "reference_expression": {
                "path":
                    str(
                        REFERENCE_EXPRESSION_FILE
                        .relative_to(
                            SOURCE_ROOT
                        )
                    ),
                "sha256":
                    sha256_file(
                        REFERENCE_EXPRESSION_FILE
                    ),
            },
            "strict_weights": {
                "path":
                    str(
                        STRICT_WEIGHTS_FILE
                        .relative_to(
                            SOURCE_ROOT
                        )
                    ),
                "sha256":
                    sha256_file(
                        STRICT_WEIGHTS_FILE
                    ),
            },
            "gene_map": {
                "path":
                    str(
                        GENE_MAP_FILE
                        .relative_to(
                            SOURCE_ROOT
                        )
                    ),
                "sha256":
                    sha256_file(
                        GENE_MAP_FILE
                    ),
            },
        },
        "fixture_contract": {
            "primary_modules":
                PRIMARY_MODULES,
            "minimum_module_genes":
                3,
            "reference_gene_standardization":
                "within_cohort_sample_sd_ddof_1",
            "external_gene_standardization":
                "within_cohort_sample_sd_ddof_1",
            "edge_matrix":
                "pearson_gene_gene_correlation",
            "edge_summary":
                (
                    "spearman_between_upper_"
                    "triangles"
                ),
            "external_loading_model":
                "first_principal_component",
            "external_pc1_orientation":
                (
                    "orient_to_positive_"
                    "correlation_with_frozen_"
                    "signed_score"
                ),
            "loading_summary":
                (
                    "spearman_frozen_risk_"
                    "loadings_vs_oriented_"
                    "external_pc1_loadings"
                ),
            "outcome_loaded":
                False,
            "initial_reusable_scope":
                (
                    "direct_preservation_"
                    "metrics_before_random_"
                    "controls"
                ),
        },
        "fixture_dimensions": {
            "reference_samples":
                int(
                    reference_fixture.shape[0]
                ),
            "external_samples":
                int(
                    external_fixture.shape[0]
                ),
            "common_fixture_genes":
                int(
                    len(genes)
                ),
            "primary_weight_rows":
                int(
                    primary_weights.shape[0]
                ),
        },
        "module_coverage":
            module_counts,
        "fixture_files":
            fixture_files,
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
        "Preservation fixture summary"
    )
    print("=" * 80)
    print(
        "Reference samples: "
        f"{reference_fixture.shape[0]}"
    )
    print(
        "External samples: "
        f"{external_fixture.shape[0]}"
    )
    print(
        "Common fixture genes: "
        f"{len(genes)}"
    )
    print(
        "Primary frozen weight rows: "
        f"{primary_weights.shape[0]}"
    )

    print("")
    print(
        coverage[
            [
                "module_label",
                "n_frozen_genes",
                "n_common_genes",
                "coverage_fraction",
            ]
        ].to_string(
            index=False
        )
    )

    print("")
    print(
        "Expected direct preservation:"
    )
    print("")

    direct_columns = [
        column
        for column in [
            "module_label",
            "edge_spearman",
            "loading_spearman",
            (
                "external_pc1_variance_"
                "explained"
            ),
            (
                "pc1_orientation_correlation_"
                "with_frozen_score"
            ),
        ]
        if column
        in structure.columns
    ]

    print(
        structure[
            direct_columns
        ].to_string(
            index=False
        )
    )

    print("")
    print(
        "Fixture manifest:"
    )
    print(
        OUTPUT_MANIFEST.relative_to(
            ROOT
        )
    )
    print("=" * 80)


if __name__ == "__main__":
    main()
