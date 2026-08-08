"""Audit exact replayability of the legacy GSE239948 random-panel control."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any
import types

import numpy as np
import pandas as pd


TOOL_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]

CORE_DIR = ROOT / "core"

if str(CORE_DIR) not in sys.path:
    sys.path.insert(
        0,
        str(CORE_DIR),
    )

from transport_audit.random_controls import (
    LEGACY_Z_VARIANCE,
    variability_bins as core_variability_bins,
)


SOURCE_ROOT = (
    ROOT.parent
    / "paper4_sarcoma_dog"
)

LEGACY_SCRIPT = (
    SOURCE_ROOT
    / "scripts"
    / "46_gse239948_external_canine_representation_v2.py"
)

LOCKED_RANDOM_FILE = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE239948_external_random_panel_controls.csv"
)

LOCKED_STRUCTURE_FILE = (
    SOURCE_ROOT
    / "results"
    / "tables"
    / "GSE239948_external_module_structure_preservation.csv"
)

SERIALIZED_EXTERNAL_FILE = (
    SOURCE_ROOT
    / "data"
    / "processed"
    / "canine_validation_GSE239948_expression_log2_symbol.csv"
)

FIXTURE_DIR = (
    ROOT
    / "reference_results"
    / "osteosarcoma_locked"
    / "preservation_fixture"
)

FIXTURE_REFERENCE = (
    FIXTURE_DIR
    / "random_panel_DOG2_prestandardization.npy"
)

FIXTURE_EXTERNAL = (
    FIXTURE_DIR
    / "random_panel_GSE239948_prestandardization.npy"
)

FIXTURE_GENES = (
    FIXTURE_DIR
    / "random_panel_background_genes.csv"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "random_panel_legacy_replay_audit"
)

OUTPUT_REPLAY = (
    OUTPUT_DIR
    / "authoritative_v2_current_environment_replay.csv"
)

OUTPUT_SUMMARY = (
    OUTPUT_DIR
    / "random_panel_legacy_replay_audit.json"
)


MODULES = [
    "M34",
    "M11",
    "M24",
    "M40",
]

RANDOM_COLUMNS = [
    "n_module_genes",
    "n_random_panels",
    "observed_edge_spearman",
    "random_edge_median",
    "random_edge_q05",
    "random_edge_q95",
    "random_panel_empirical_p",
]


def require_path(
    path: Path,
) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"Required file not found:\n{path}"
        )


def load_legacy_module():
    """
    Load authoritative script 46 v2 without requiring its plotting stack.

    The audit uses only scientific/data-processing functions from the
    legacy module. matplotlib is imported by the legacy script solely
    because it also contains figure-generation code.
    """
    specification = (
        importlib.util.spec_from_file_location(
            "legacy_script46_v2",
            LEGACY_SCRIPT,
        )
    )

    if (
        specification is None
        or specification.loader is None
    ):
        raise RuntimeError(
            "Could not load legacy script 46 v2."
        )

    injected_modules: list[str] = []

    if "matplotlib" not in sys.modules:
        matplotlib_stub = types.ModuleType(
            "matplotlib"
        )

        matplotlib_stub.__path__ = []

        pyplot_stub = types.ModuleType(
            "matplotlib.pyplot"
        )

        matplotlib_stub.pyplot = (
            pyplot_stub
        )

        sys.modules[
            "matplotlib"
        ] = matplotlib_stub

        sys.modules[
            "matplotlib.pyplot"
        ] = pyplot_stub

        injected_modules.extend(
            [
                "matplotlib.pyplot",
                "matplotlib",
            ]
        )

    module = (
        importlib.util.module_from_spec(
            specification
        )
    )

    sys.modules[
        specification.name
    ] = module

    try:
        specification.loader.exec_module(
            module
        )
    except Exception:
        sys.modules.pop(
            specification.name,
            None,
        )
        raise
    finally:
        for name in injected_modules:
            sys.modules.pop(
                name,
                None,
            )

    return module


def find_cached_raw_file() -> Path:
    filename = (
        "GSE239948_CCOGC.txt.gz"
    )

    candidates = [
        (
            SOURCE_ROOT
            / "data"
            / "raw"
            / "canine_validation_GSE239948"
            / filename
        ),
        (
            SOURCE_ROOT
            / "data"
            / "raw"
            / "GSE239948"
            / filename
        ),
        (
            SOURCE_ROOT
            / "data"
            / "raw"
            / filename
        ),
        SOURCE_ROOT / filename,
    ]

    for path in candidates:
        if path.exists():
            return path

    recursive = list(
        SOURCE_ROOT.rglob(
            filename
        )
    )

    if recursive:
        return recursive[0]

    raise FileNotFoundError(
        "Cached raw GSE239948 file "
        "was not found. No download "
        "was attempted."
    )


def order_modules(
    table: pd.DataFrame,
) -> pd.DataFrame:
    result = table.copy()

    result[
        "module_label"
    ] = (
        result[
            "module_label"
        ]
        .astype(str)
    )

    return (
        result.set_index(
            "module_label"
        )
        .reindex(
            MODULES
        )
        .reset_index()
    )


def compare_random_tables(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
) -> tuple[
    bool,
    float,
    list[str],
]:
    observed = order_modules(
        observed
    )

    expected = order_modules(
        expected
    )

    mismatches: list[str] = []
    maximum_difference = 0.0

    for column in RANDOM_COLUMNS:
        left = pd.to_numeric(
            observed[column],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        right = pd.to_numeric(
            expected[column],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        differences = np.abs(
            left - right
        )

        if np.isfinite(
            differences
        ).any():
            maximum_difference = max(
                maximum_difference,
                float(
                    np.nanmax(
                        differences
                    )
                ),
            )

        if not np.allclose(
            left,
            right,
            rtol=1.0e-8,
            atol=1.0e-10,
            equal_nan=True,
        ):
            mismatches.append(
                column
            )

    return (
        len(mismatches) == 0,
        maximum_difference,
        mismatches,
    )


def bin_agreement(
    left: pd.Series,
    right: pd.Series,
) -> dict[str, Any]:
    common = (
        left.index
        .intersection(
            right.index
        )
    )

    left_values = (
        left.reindex(
            common
        ).to_numpy()
    )

    right_values = (
        right.reindex(
            common
        ).to_numpy()
    )

    both_nan = (
        pd.isna(
            left_values
        )
        & pd.isna(
            right_values
        )
    )

    same = (
        left_values
        == right_values
    ) | both_nan

    return {
        "n_common":
            int(
                len(common)
            ),
        "exact_agreement_fraction":
            float(
                np.mean(
                    same
                )
            ),
    }


def matrix_difference(
    left: pd.DataFrame,
    right: pd.DataFrame,
) -> dict[str, Any]:
    common_columns = (
        left.columns
        .intersection(
            right.columns
        )
    )

    common_rows = (
        left.index
        .intersection(
            right.index
        )
    )

    left_values = (
        left.loc[
            common_rows,
            common_columns,
        ].to_numpy(
            dtype=float
        )
    )

    right_values = (
        right.loc[
            common_rows,
            common_columns,
        ].to_numpy(
            dtype=float
        )
    )

    differences = np.abs(
        left_values
        - right_values
    )

    return {
        "rows":
            int(
                len(
                    common_rows
                )
            ),
        "columns":
            int(
                len(
                    common_columns
                )
            ),
        "maximum_absolute_difference":
            float(
                np.nanmax(
                    differences
                )
            ),
        "mean_absolute_difference":
            float(
                np.nanmean(
                    differences
                )
            ),
    }


def reconstruct_authoritative_inputs(
    legacy,
) -> dict[str, Any]:
    raw_path = (
        find_cached_raw_file()
    )

    reference = pd.read_csv(
        legacy.REFERENCE_EXPRESSION_FILE,
        index_col=0,
        low_memory=False,
    )

    weights = pd.read_csv(
        legacy.STRICT_WEIGHTS_FILE,
        low_memory=False,
    )

    raw_table = legacy.read_raw_table(
        raw_path
    )

    samples = legacy.sample_columns(
        raw_table
    )

    identifier_audit = (
        legacy.identifier_pair_audit(
            table=raw_table,
            sample_cols=samples,
            weights=weights,
        )
    )

    (
        gene_column,
        frozen_identifier_column,
        identifier_scheme,
    ) = legacy.choose_identifier_pair(
        identifier_audit
    )

    weights = weights.copy()

    weights[
        "canine_gene_symbol"
    ] = weights[
        frozen_identifier_column
    ].map(
        lambda value:
            legacy.normalize_identifier(
                value,
                identifier_scheme,
            )
    )

    (
        external_raw,
        _,
    ) = legacy.collapse_to_identifiers(
        table=raw_table,
        gene_column=gene_column,
        sample_cols=samples,
        identifier_scheme=(
            identifier_scheme
        ),
    )

    (
        external,
        transform_diagnostics,
    ) = legacy.choose_transform(
        external_raw
    )

    (
        reference,
        _,
    ) = (
        legacy.remap_reference_expression(
            reference=reference,
            weights=weights,
            frozen_identifier_column=(
                frozen_identifier_column
            ),
            identifier_scheme=(
                identifier_scheme
            ),
        )
    )

    external = external.loc[
        :,
        ~pd.Index(
            external.columns
        ).duplicated(
            keep="first"
        ),
    ].copy()

    reference_z = (
        legacy.zscore_columns(
            reference
        )
    )

    external_z = (
        legacy.zscore_columns(
            external
        )
    )

    return {
        "raw_path":
            raw_path,
        "reference":
            reference,
        "external":
            external,
        "reference_z":
            reference_z,
        "external_z":
            external_z,
        "weights":
            weights,
        "gene_column":
            gene_column,
        "frozen_identifier_column":
            frozen_identifier_column,
        "identifier_scheme":
            identifier_scheme,
        "transform_diagnostics":
            transform_diagnostics,
    }


def main() -> None:
    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- legacy random-panel replay audit"
    )
    print("=" * 80)
    print(
        f"Tool version: "
        f"{TOOL_VERSION}"
    )

    for path in [
        LEGACY_SCRIPT,
        LOCKED_RANDOM_FILE,
        LOCKED_STRUCTURE_FILE,
        SERIALIZED_EXTERNAL_FILE,
        FIXTURE_REFERENCE,
        FIXTURE_EXTERNAL,
        FIXTURE_GENES,
    ]:
        require_path(
            path
        )

    legacy = load_legacy_module()

    if (
        legacy.SCRIPT_VERSION
        != (
            "46-gse239948-"
            "external-canine-"
            "representation-v2"
        )
    ):
        raise RuntimeError(
            "Loaded script is not "
            "authoritative script 46 v2."
        )

    print("")
    print(
        "Reconstructing the exact "
        "authoritative script-46 inputs..."
    )

    inputs = (
        reconstruct_authoritative_inputs(
            legacy
        )
    )

    print(
        "Selected mapping:"
    )
    print(
        "  External column: "
        f"{inputs['gene_column']}"
    )
    print(
        "  Frozen identifier column: "
        f"{inputs['frozen_identifier_column']}"
    )
    print(
        "  Identifier scheme: "
        f"{inputs['identifier_scheme']}"
    )
    print(
        "  Transform: "
        f"{inputs['transform_diagnostics']['transform']}"
    )

    reference_z = inputs[
        "reference_z"
    ]

    external_z = inputs[
        "external_z"
    ]

    legacy_common = (
        reference_z.columns
        .intersection(
            external_z.columns
        )
    )

    print("")
    print(
        "Authoritative in-memory universe:"
    )
    print(
        "  Reference samples: "
        f"{reference_z.shape[0]}"
    )
    print(
        "  External samples: "
        f"{external_z.shape[0]}"
    )
    print(
        "  Shared genes: "
        f"{len(legacy_common)}"
    )

    locked_structure = pd.read_csv(
        LOCKED_STRUCTURE_FILE
    )

    locked_random = pd.read_csv(
        LOCKED_RANDOM_FILE
    )

    print("")
    print("=" * 80)
    print(
        "Replay authoritative v2 "
        "random_panel_controls()"
    )
    print("=" * 80)
    print(
        "This executes the original "
        "script-46 function unchanged."
    )
    print("")

    replay = legacy.random_panel_controls(
        reference_z=reference_z,
        external_z=external_z,
        weights=inputs[
            "weights"
        ],
        observed=locked_structure,
    )

    (
        replay_matches_lock,
        replay_max_difference,
        replay_mismatches,
    ) = compare_random_tables(
        replay,
        locked_random,
    )

    print(
        "Authoritative current-environment "
        "replay vs locked output: "
        + (
            "PASS"
            if replay_matches_lock
            else "FAIL"
        )
    )

    print(
        "Maximum absolute difference: "
        f"{replay_max_difference:.6e}"
    )

    if replay_mismatches:
        print(
            "Mismatched columns: "
            + ", ".join(
                replay_mismatches
            )
        )

    print("")
    print(
        replay.to_string(
            index=False
        )
    )

    raw_legacy_bins = (
        legacy.variability_bins(
            reference_z,
            external_z,
        )
    )

    serialized_external = (
        pd.read_csv(
            SERIALIZED_EXTERNAL_FILE,
            index_col=0,
            low_memory=False,
        )
    )

    serialized_external.index = (
        serialized_external.index
        .astype(str)
    )

    serialized_external_z = (
        legacy.zscore_columns(
            serialized_external
        )
    )

    serialized_bins = (
        legacy.variability_bins(
            reference_z,
            serialized_external_z,
        )
    )

    genes = (
        pd.read_csv(
            FIXTURE_GENES
        )[
            "gene_symbol"
        ]
        .astype(str)
        .tolist()
    )

    fixture_reference = (
        pd.DataFrame(
            np.load(
                FIXTURE_REFERENCE,
                allow_pickle=False,
            ),
            columns=genes,
        )
    )

    fixture_external = (
        pd.DataFrame(
            np.load(
                FIXTURE_EXTERNAL,
                allow_pickle=False,
            ),
            columns=genes,
        )
    )

    fixture_reference_z = (
        legacy.zscore_columns(
            fixture_reference
        )
    )

    fixture_external_z = (
        legacy.zscore_columns(
            fixture_external
        )
    )

    fixture_legacy_bins = (
        legacy.variability_bins(
            fixture_reference_z,
            fixture_external_z,
        )
    )

    core_raw_bins = (
        core_variability_bins(
            inputs[
                "reference"
            ],
            inputs[
                "external"
            ],
            matching_mode=(
                LEGACY_Z_VARIANCE
            ),
            n_bins=10,
        )
    )

    raw_vs_serialized_bins = (
        bin_agreement(
            raw_legacy_bins,
            serialized_bins,
        )
    )

    raw_vs_fixture_bins = (
        bin_agreement(
            raw_legacy_bins,
            fixture_legacy_bins,
        )
    )

    raw_vs_core_bins = (
        bin_agreement(
            raw_legacy_bins,
            core_raw_bins,
        )
    )

    external_z_difference = (
        matrix_difference(
            external_z,
            serialized_external_z,
        )
    )

    legacy_common_list = [
        str(gene)
        for gene in legacy_common
    ]

    fixture_gene_set_equal = bool(
        set(
            legacy_common_list
        )
        == set(
            genes
        )
    )

    fixture_gene_order_equal = bool(
        legacy_common_list
        == genes
    )

    print("")
    print("=" * 80)
    print(
        "Numerical-path audit"
    )
    print("=" * 80)

    print(
        "Legacy shared gene set equals "
        "fixture gene set: "
        f"{fixture_gene_set_equal}"
    )

    print(
        "Legacy shared gene order equals "
        "fixture gene order: "
        f"{fixture_gene_order_equal}"
    )

    print("")
    print(
        "In-memory external z-score vs "
        "serialized-output z-score:"
    )

    print(
        "  Maximum absolute difference: "
        f"{external_z_difference['maximum_absolute_difference']:.3e}"
    )

    print(
        "  Mean absolute difference: "
        f"{external_z_difference['mean_absolute_difference']:.3e}"
    )

    print("")
    print(
        "Exact bin agreement:"
    )

    print(
        "  Authoritative in-memory vs "
        "serialized external: "
        f"{raw_vs_serialized_bins['exact_agreement_fraction']:.3%}"
    )

    print(
        "  Authoritative in-memory vs "
        "fixture replay: "
        f"{raw_vs_fixture_bins['exact_agreement_fraction']:.3%}"
    )

    print(
        "  Authoritative legacy function vs "
        "core on same raw inputs: "
        f"{raw_vs_core_bins['exact_agreement_fraction']:.3%}"
    )

    if not replay_matches_lock:
        diagnosis = (
            "authoritative_function_not_"
            "reproducible_against_locked_"
            "random_output_in_current_"
            "environment"
        )
    elif (
        raw_vs_core_bins[
            "exact_agreement_fraction"
        ]
        < 1.0
    ):
        diagnosis = (
            "core_legacy_bin_"
            "implementation_difference"
        )
    elif (
        raw_vs_fixture_bins[
            "exact_agreement_fraction"
        ]
        < 1.0
    ):
        diagnosis = (
            "fixture_or_serialization_"
            "changes_roundoff_defined_"
            "legacy_bins"
        )
    else:
        diagnosis = (
            "sampling_path_requires_"
            "further_audit"
        )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    replay.to_csv(
        OUTPUT_REPLAY,
        index=False,
    )

    summary = {
        "tool_version":
            TOOL_VERSION,
        "authoritative_script_version":
            legacy.SCRIPT_VERSION,
        "selected_mapping": {
            "external_column":
                inputs[
                    "gene_column"
                ],
            "frozen_identifier_column":
                inputs[
                    "frozen_identifier_column"
                ],
            "identifier_scheme":
                inputs[
                    "identifier_scheme"
                ],
            "transform":
                inputs[
                    "transform_diagnostics"
                ][
                    "transform"
                ],
        },
        "dimensions": {
            "reference_samples":
                int(
                    reference_z.shape[0]
                ),
            "external_samples":
                int(
                    external_z.shape[0]
                ),
            "shared_genes":
                int(
                    len(
                        legacy_common
                    )
                ),
        },
        "authoritative_replay": {
            "matches_locked_output":
                replay_matches_lock,
            "maximum_absolute_difference":
                replay_max_difference,
            "mismatched_columns":
                replay_mismatches,
        },
        "gene_universe": {
            "fixture_set_equal":
                fixture_gene_set_equal,
            "fixture_order_equal":
                fixture_gene_order_equal,
        },
        "external_serialization": {
            "z_matrix_difference":
                external_z_difference,
        },
        "bin_agreement": {
            (
                "authoritative_in_memory_"
                "vs_serialized"
            ):
                raw_vs_serialized_bins[
                    "exact_agreement_fraction"
                ],
            (
                "authoritative_in_memory_"
                "vs_fixture"
            ):
                raw_vs_fixture_bins[
                    "exact_agreement_fraction"
                ],
            (
                "authoritative_legacy_"
                "vs_core_same_raw_inputs"
            ):
                raw_vs_core_bins[
                    "exact_agreement_fraction"
                ],
        },
        "diagnosis":
            diagnosis,
    }

    OUTPUT_SUMMARY.write_text(
        json.dumps(
            summary,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("")
    print("=" * 80)
    print(
        "Audit diagnosis"
    )
    print("=" * 80)

    print(
        diagnosis
    )

    print("")
    print(
        "Summary: "
        + str(
            OUTPUT_SUMMARY
            .relative_to(ROOT)
        )
    )

    print("=" * 80)


if __name__ == "__main__":
    main()
