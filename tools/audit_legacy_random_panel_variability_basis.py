"""Audit the variability basis of the legacy GSE239948 random-panel control."""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats


TOOL_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]

SOURCE_ROOT = (
    ROOT.parent
    / "paper4_sarcoma_dog"
)

REFERENCE_FILE = (
    SOURCE_ROOT
    / "data"
    / "processed"
    / "GSE238110_DOG2_expression_log2cpm_matched_allgenes.csv"
)

EXTERNAL_FILE = (
    SOURCE_ROOT
    / "data"
    / "processed"
    / "canine_validation_GSE239948_expression_log2_symbol.csv"
)

GENE_MAP_FILE = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE239948_external_gene_mapping.csv"
)

WEIGHTS_FILE = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE238110_frozen_transfer_gene_weights_strict.csv"
)

LEGACY_RANDOM_FILE = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE239948_external_random_panel_controls.csv"
)

LEGACY_CLASSIFICATION_FILE = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE239948_external_representation_classification.csv"
)

LEGACY_MANIFEST_FILE = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE239948_external_representation_manifest.json"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "random_panel_matching_audit"
)

GENE_AUDIT_FILE = (
    OUTPUT_DIR
    / "gene_variability_and_bin_audit.csv"
)

MODULE_AUDIT_FILE = (
    OUTPUT_DIR
    / "module_variability_bin_audit.csv"
)

SUMMARY_FILE = (
    OUTPUT_DIR
    / "random_panel_variability_basis_audit.json"
)


PRIMARY_MODULES = [
    "M34",
    "M11",
    "M24",
    "M40",
]

N_VARIABILITY_BINS = 10

EXPECTED_COMMON_MODULE_COUNTS = {
    "M34": 154,
    "M11": 6,
    "M24": 6,
    "M40": 106,
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


def require_path(
    path: Path,
) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"Required source file not found:\n"
            f"{path}"
        )


def clean_symbol(
    value: Any,
) -> str:
    text = str(value).strip().upper()

    if text in {
        "",
        "NAN",
        "NONE",
        "NA",
    }:
        return ""

    if "|" in text:
        parts = [
            part.strip()
            for part
            in text.split("|")
            if part.strip()
        ]

        symbol_like = [
            part
            for part in parts
            if not part.startswith(
                "ENSCAFG"
            )
        ]

        if parts:
            text = (
                symbol_like[-1]
                if symbol_like
                else parts[-1]
            )

    text = re.sub(
        r"\.\d+$",
        "",
        text,
    )

    return text


def normalize_identifier_symbol(
    value: Any,
) -> str:
    """Mirror the symbol normalization used by script 46 v2."""
    text = str(value).strip().upper()

    if text in {
        "",
        "NAN",
        "NONE",
        "NA",
        "---",
    }:
        return ""

    gene_name_match = re.search(
        (
            r"(?:GENE[_ ]?NAME|SYMBOL)"
            r"\s*[=:]\s*"
            r"([A-Z0-9_.-]+)"
        ),
        text,
    )

    if gene_name_match:
        text = gene_name_match.group(1)

    if "|" in text:
        parts = [
            part.strip()
            for part
            in text.split("|")
            if part.strip()
        ]

        non_ensembl = [
            part
            for part in parts
            if not re.fullmatch(
                r"ENSCAFG\d+(?:\.\d+)?",
                part,
            )
        ]

        text = (
            non_ensembl[-1]
            if non_ensembl
            else ""
        )

    if ";" in text:
        parts = [
            part.strip()
            for part
            in text.split(";")
            if part.strip()
        ]

        symbol_like = [
            part
            for part in parts
            if (
                re.fullmatch(
                    r"[A-Z][A-Z0-9_.-]*",
                    part,
                )
                and not part.startswith(
                    "ENSCAFG"
                )
            )
        ]

        if symbol_like:
            text = symbol_like[-1]

    text = re.sub(
        r"\.\d+$",
        "",
        text,
    )

    text = re.sub(
        r"_\d+$",
        "",
        text,
    )

    text = text.strip()

    if text.startswith(
        "ENSCAFG"
    ):
        return ""

    if not re.fullmatch(
        r"[A-Z][A-Z0-9_.-]*",
        text,
    ):
        return ""

    return text


def reconstruct_reference(
    reference: pd.DataFrame,
    gene_map: pd.DataFrame,
) -> pd.DataFrame:
    required = {
        "analysis_gene_id",
        "reference_expression_column",
    }

    missing = sorted(
        required.difference(
            gene_map.columns
        )
    )

    if missing:
        raise ValueError(
            "Gene map is missing columns: "
            + ", ".join(missing)
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

    mapping = mapping.drop_duplicates(
        "reference_expression_column",
        keep="first",
    )

    exact_map = dict(
        zip(
            mapping[
                "reference_expression_column"
            ],
            mapping[
                "analysis_gene_id"
            ],
        )
    )

    rename_map: dict[
        str,
        str,
    ] = {}

    for column in (
        reference.columns
        .astype(str)
    ):
        analysis_id = (
            exact_map.get(
                column,
                "",
            )
        )

        if not analysis_id:
            analysis_id = (
                normalize_identifier_symbol(
                    column
                )
            )

        if not analysis_id:
            continue

        rename_map[
            column
        ] = analysis_id

    remapped = reference.rename(
        columns=rename_map
    )

    remapped = remapped.loc[
        :,
        list(
            rename_map.values()
        ),
    ]

    remapped = remapped.loc[
        :,
        ~pd.Index(
            remapped.columns
        ).duplicated(
            keep="first"
        ),
    ].copy()

    return remapped


def prepare_numeric(
    expression: pd.DataFrame,
) -> pd.DataFrame:
    x = expression.apply(
        pd.to_numeric,
        errors="coerce",
    )

    x = x.fillna(
        x.median(axis=0)
    )

    return x


def zscore_columns(
    expression: pd.DataFrame,
) -> pd.DataFrame:
    x = prepare_numeric(
        expression
    )

    std = (
        x.std(
            axis=0,
            ddof=1,
        )
        .replace(
            0.0,
            np.nan,
        )
    )

    z = (
        x
        - x.mean(axis=0)
    ) / std

    return z.loc[
        :,
        z.notna().all(
            axis=0
        ),
    ].copy()


def variability_bins(
    reference_expression: pd.DataFrame,
    external_expression: pd.DataFrame,
) -> tuple[
    pd.Series,
    pd.Series,
]:
    common = (
        reference_expression.columns
        .intersection(
            external_expression.columns
        )
    )

    reference_variance = (
        reference_expression[
            common
        ].var(
            axis=0,
            ddof=1,
        )
    )

    external_variance = (
        external_expression[
            common
        ].var(
            axis=0,
            ddof=1,
        )
    )

    combined_rank = (
        reference_variance.rank(
            pct=True
        )
        + external_variance.rank(
            pct=True
        )
    ) / 2.0

    bins = pd.qcut(
        combined_rank.rank(
            method="average"
        ),
        q=min(
            N_VARIABILITY_BINS,
            len(combined_rank),
        ),
        labels=False,
        duplicates="drop",
    )

    return (
        bins,
        combined_rank,
    )


def variance_summary(
    values: pd.Series,
) -> dict[str, Any]:
    numeric = pd.to_numeric(
        values,
        errors="coerce",
    ).dropna()

    absolute_deviation = (
        numeric
        - 1.0
    ).abs()

    return {
        "n":
            int(
                numeric.shape[0]
            ),
        "minimum":
            float(
                numeric.min()
            ),
        "median":
            float(
                numeric.median()
            ),
        "maximum":
            float(
                numeric.max()
            ),
        "standard_deviation":
            float(
                numeric.std(
                    ddof=1
                )
            ),
        "max_abs_deviation_from_one":
            float(
                absolute_deviation.max()
            ),
        "fraction_within_1e_12_of_one":
            float(
                (
                    absolute_deviation
                    <= 1.0e-12
                ).mean()
            ),
        "fraction_within_1e_10_of_one":
            float(
                (
                    absolute_deviation
                    <= 1.0e-10
                ).mean()
            ),
        "unique_exact":
            int(
                numeric.nunique(
                    dropna=True
                )
            ),
        "unique_rounded_12dp":
            int(
                numeric.round(
                    12
                ).nunique(
                    dropna=True
                )
            ),
    }


def main() -> None:
    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- legacy random-panel variability audit"
    )
    print("=" * 80)

    for path in [
        REFERENCE_FILE,
        EXTERNAL_FILE,
        GENE_MAP_FILE,
        WEIGHTS_FILE,
        LEGACY_RANDOM_FILE,
        LEGACY_CLASSIFICATION_FILE,
        LEGACY_MANIFEST_FILE,
    ]:
        require_path(path)

    legacy_manifest = json.loads(
        LEGACY_MANIFEST_FILE.read_text(
            encoding="utf-8"
        )
    )

    if (
        str(
            legacy_manifest.get(
                "script_version",
                "",
            )
        )
        != (
            "46-gse239948-"
            "external-canine-"
            "representation-v2"
        )
    ):
        raise RuntimeError(
            "Unexpected legacy script "
            "version."
        )

    if bool(
        legacy_manifest.get(
            "outcome_loaded",
            True,
        )
    ):
        raise RuntimeError(
            "Legacy manifest does not "
            "confirm outcome_loaded=false."
        )

    reference_raw = pd.read_csv(
        REFERENCE_FILE,
        index_col=0,
        low_memory=False,
    )

    external_raw = pd.read_csv(
        EXTERNAL_FILE,
        index_col=0,
        low_memory=False,
    )

    gene_map = pd.read_csv(
        GENE_MAP_FILE,
        low_memory=False,
    )

    weights = pd.read_csv(
        WEIGHTS_FILE,
        low_memory=False,
    )

    legacy_random = pd.read_csv(
        LEGACY_RANDOM_FILE
    )

    legacy_classification = pd.read_csv(
        LEGACY_CLASSIFICATION_FILE
    )

    reference = reconstruct_reference(
        reference_raw,
        gene_map,
    )

    external = external_raw.copy()

    external.columns = (
        external.columns
        .astype(str)
        .str.strip()
    )

    external = external.loc[
        :,
        ~pd.Index(
            external.columns
        ).duplicated(
            keep="first"
        ),
    ].copy()

    reference_numeric = prepare_numeric(
        reference
    )

    external_numeric = prepare_numeric(
        external
    )

    reference_z = zscore_columns(
        reference_numeric
    )

    external_z = zscore_columns(
        external_numeric
    )

    common = (
        reference_z.columns
        .intersection(
            external_z.columns
        )
    )

    if len(common) < 100:
        raise RuntimeError(
            "Unexpectedly small shared "
            "background universe."
        )

    reference_raw_common = (
        reference_numeric[
            common
        ]
    )

    external_raw_common = (
        external_numeric[
            common
        ]
    )

    reference_z_common = (
        reference_z[
            common
        ]
    )

    external_z_common = (
        external_z[
            common
        ]
    )

    reference_raw_variance = (
        reference_raw_common.var(
            axis=0,
            ddof=1,
        )
    )

    external_raw_variance = (
        external_raw_common.var(
            axis=0,
            ddof=1,
        )
    )

    reference_z_variance = (
        reference_z_common.var(
            axis=0,
            ddof=1,
        )
    )

    external_z_variance = (
        external_z_common.var(
            axis=0,
            ddof=1,
        )
    )

    (
        legacy_bins,
        legacy_combined_rank,
    ) = variability_bins(
        reference_z_common,
        external_z_common,
    )

    (
        prestandardization_bins,
        prestandardization_combined_rank,
    ) = variability_bins(
        reference_raw_common,
        external_raw_common,
    )

    weights = weights.copy()

    weights[
        "canine_gene_symbol"
    ] = (
        weights[
            "canine_gene_symbol"
        ]
        .map(
            clean_symbol
        )
    )

    frozen_union = set(
        weights.loc[
            weights[
                "module_label"
            ].isin(
                PRIMARY_MODULES
            ),
            "canine_gene_symbol",
        ]
    )

    frozen_union.discard("")

    candidate_genes = [
        gene
        for gene in common
        if gene
        not in frozen_union
    ]

    gene_audit = pd.DataFrame(
        {
            "gene":
                common,
            "reference_raw_variance":
                reference_raw_variance
                .reindex(
                    common
                )
                .to_numpy(),
            "external_raw_variance":
                external_raw_variance
                .reindex(
                    common
                )
                .to_numpy(),
            "reference_z_variance":
                reference_z_variance
                .reindex(
                    common
                )
                .to_numpy(),
            "external_z_variance":
                external_z_variance
                .reindex(
                    common
                )
                .to_numpy(),
            "legacy_combined_rank":
                legacy_combined_rank
                .reindex(
                    common
                )
                .to_numpy(),
            "prestandardization_combined_rank":
                prestandardization_combined_rank
                .reindex(
                    common
                )
                .to_numpy(),
            "legacy_variability_bin":
                legacy_bins
                .reindex(
                    common
                )
                .to_numpy(),
            "prestandardization_variability_bin":
                prestandardization_bins
                .reindex(
                    common
                )
                .to_numpy(),
            "is_primary_frozen_gene": [
                gene in frozen_union
                for gene in common
            ],
            "is_random_panel_candidate": [
                gene
                in candidate_genes
                for gene in common
            ],
        }
    )

    module_rows: list[
        dict[str, Any]
    ] = []

    for module in PRIMARY_MODULES:
        part = weights[
            weights[
                "module_label"
            ]
            .astype(str)
            .eq(module)
        ].copy()

        part = part.drop_duplicates(
            "canine_gene_symbol",
            keep="first",
        )

        genes = [
            gene
            for gene
            in part[
                "canine_gene_symbol"
            ]
            if gene in common
        ]

        expected_count = (
            EXPECTED_COMMON_MODULE_COUNTS[
                module
            ]
        )

        if (
            len(genes)
            != expected_count
        ):
            raise RuntimeError(
                "Shared-gene reconstruction "
                "does not match locked "
                f"coverage for {module}: "
                f"expected={expected_count}, "
                f"observed={len(genes)}"
            )

        legacy_target_bins = (
            legacy_bins
            .reindex(
                genes
            )
        )

        corrected_target_bins = (
            prestandardization_bins
            .reindex(
                genes
            )
        )

        same = (
            legacy_target_bins
            .to_numpy()
            == corrected_target_bins
            .to_numpy()
        )

        module_rows.append(
            {
                "module_label":
                    module,
                "n_target_genes":
                    len(genes),
                "legacy_bins_used":
                    int(
                        legacy_target_bins
                        .nunique(
                            dropna=True
                        )
                    ),
                "prestandardization_bins_used":
                    int(
                        corrected_target_bins
                        .nunique(
                            dropna=True
                        )
                    ),
                "fraction_same_bin":
                    float(
                        np.mean(
                            same
                        )
                    ),
                "legacy_bin_counts":
                    json.dumps(
                        legacy_target_bins
                        .value_counts(
                            sort=False
                        )
                        .sort_index()
                        .astype(int)
                        .to_dict(),
                        sort_keys=True,
                    ),
                "prestandardization_bin_counts":
                    json.dumps(
                        corrected_target_bins
                        .value_counts(
                            sort=False
                        )
                        .sort_index()
                        .astype(int)
                        .to_dict(),
                        sort_keys=True,
                    ),
            }
        )

    module_audit = pd.DataFrame(
        module_rows
    )

    rank_correlation = (
        stats.spearmanr(
            legacy_combined_rank
            .reindex(common)
            .to_numpy(
                dtype=float
            ),
            prestandardization_combined_rank
            .reindex(common)
            .to_numpy(
                dtype=float
            ),
        ).statistic
    )

    bin_agreement = float(
        np.mean(
            legacy_bins
            .reindex(common)
            .to_numpy()
            == prestandardization_bins
            .reindex(common)
            .to_numpy()
        )
    )

    reference_z_summary = (
        variance_summary(
            reference_z_variance
        )
    )

    external_z_summary = (
        variance_summary(
            external_z_variance
        )
    )

    reference_raw_summary = (
        variance_summary(
            reference_raw_variance
        )
    )

    external_raw_summary = (
        variance_summary(
            external_raw_variance
        )
    )

    z_variance_effectively_constant = bool(
        (
            reference_z_variance
            .sub(1.0)
            .abs()
            .max()
            <= 1.0e-10
        )
        and (
            external_z_variance
            .sub(1.0)
            .abs()
            .max()
            <= 1.0e-10
        )
    )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    gene_audit.to_csv(
        GENE_AUDIT_FILE,
        index=False,
    )

    module_audit.to_csv(
        MODULE_AUDIT_FILE,
        index=False,
    )

    summary = {
        "tool_version":
            TOOL_VERSION,
        "generated_at_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),
        "source_script_version":
            legacy_manifest[
                "script_version"
            ],
        "outcome_loaded":
            False,
        "n_reference_samples":
            int(
                reference_z.shape[0]
            ),
        "n_external_samples":
            int(
                external_z.shape[0]
            ),
        "n_shared_genes":
            int(
                len(common)
            ),
        "n_primary_frozen_genes_in_shared_universe":
            int(
                sum(
                    gene in frozen_union
                    for gene in common
                )
            ),
        "n_random_panel_candidate_genes":
            int(
                len(
                    candidate_genes
                )
            ),
        "legacy_matching_basis":
            (
                "variance_of_within_cohort_"
                "z_standardized_expression"
            ),
        "reference_z_variance":
            reference_z_summary,
        "external_z_variance":
            external_z_summary,
        "reference_prestandardization_variance":
            reference_raw_summary,
        "external_prestandardization_variance":
            external_raw_summary,
        "z_variance_effectively_constant":
            z_variance_effectively_constant,
        (
            "spearman_legacy_vs_"
            "prestandardization_"
            "combined_variability_rank"
        ):
            (
                float(
                    rank_correlation
                )
                if np.isfinite(
                    rank_correlation
                )
                else None
            ),
        "overall_bin_agreement_fraction":
            bin_agreement,
        "legacy_random_panel_results":
            legacy_random.to_dict(
                orient="records"
            ),
        "legacy_classification": (
            legacy_classification[
                [
                    "module_label",
                    (
                        "random_panel_"
                        "empirical_p"
                    ),
                    (
                        "external_canine_"
                        "representation_class"
                    ),
                ]
            ].to_dict(
                orient="records"
            )
        ),
        "interpretation_guardrail": (
            "This audit does not change any "
            "frozen program, outcome, direct "
            "preservation statistic, or "
            "classification. It assesses "
            "whether the legacy random-panel "
            "matching variable represented "
            "biological variability after "
            "within-cohort z-standardization."
        ),
        "files": {
            "gene_audit": {
                "path":
                    str(
                        GENE_AUDIT_FILE
                        .relative_to(
                            ROOT
                        )
                    ),
                "sha256":
                    sha256_file(
                        GENE_AUDIT_FILE
                    ),
            },
            "module_audit": {
                "path":
                    str(
                        MODULE_AUDIT_FILE
                        .relative_to(
                            ROOT
                        )
                    ),
                "sha256":
                    sha256_file(
                        MODULE_AUDIT_FILE
                    ),
            },
        },
    }

    SUMMARY_FILE.write_text(
        json.dumps(
            summary,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("")
    print(
        "Reconstructed background universe:"
    )
    print(
        f"  Reference samples: "
        f"{reference_z.shape[0]}"
    )
    print(
        f"  External samples: "
        f"{external_z.shape[0]}"
    )
    print(
        f"  Shared genes: "
        f"{len(common)}"
    )
    print(
        f"  Candidate genes after frozen "
        f"exclusion: "
        f"{len(candidate_genes)}"
    )

    print("")
    print("=" * 80)
    print(
        "Variance after legacy z-standardization"
    )
    print("=" * 80)

    variance_table = pd.DataFrame(
        [
            {
                "cohort":
                    "DOG2",
                **reference_z_summary,
            },
            {
                "cohort":
                    "GSE239948",
                **external_z_summary,
            },
        ]
    )

    print(
        variance_table.to_string(
            index=False
        )
    )

    print("")
    print(
        "Z-standardized variance "
        "effectively constant: "
        f"{z_variance_effectively_constant}"
    )

    print("")
    print("=" * 80)
    print(
        "Pre-standardization variance"
    )
    print("=" * 80)

    raw_variance_table = pd.DataFrame(
        [
            {
                "cohort":
                    "DOG2",
                **reference_raw_summary,
            },
            {
                "cohort":
                    "GSE239948",
                **external_raw_summary,
            },
        ]
    )

    print(
        raw_variance_table.to_string(
            index=False
        )
    )

    print("")
    print("=" * 80)
    print(
        "Legacy versus pre-standardization bins"
    )
    print("=" * 80)

    print(
        "Spearman correlation of combined "
        "variability ranks: "
        f"{rank_correlation:.6f}"
    )

    print(
        "Overall exact bin agreement: "
        f"{bin_agreement:.3%}"
    )

    print("")
    print(
        module_audit.to_string(
            index=False
        )
    )

    print("")
    print("=" * 80)
    print(
        "Locked legacy random-panel results"
    )
    print("=" * 80)

    print(
        legacy_random.to_string(
            index=False
        )
    )

    print("")
    print(
        "Audit summary:"
    )
    print(
        SUMMARY_FILE.relative_to(
            ROOT
        )
    )

    print("=" * 80)


if __name__ == "__main__":
    main()
