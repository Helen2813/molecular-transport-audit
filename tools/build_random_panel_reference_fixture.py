"""Build random-panel fixture from the authoritative script-46 raw path."""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from audit_random_panel_legacy_replay import (
    load_legacy_module,
    reconstruct_authoritative_inputs,
)


TOOL_VERSION = "0.2.0"

ROOT = Path(__file__).resolve().parents[1]

SOURCE_ROOT = (
    ROOT.parent
    / "paper4_sarcoma_dog"
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

FIXTURE_DIR = (
    ROOT
    / "reference_results"
    / "osteosarcoma_locked"
    / "preservation_fixture"
)

OUTPUT_REFERENCE = (
    FIXTURE_DIR
    / "random_panel_DOG2_prestandardization.npy"
)

OUTPUT_EXTERNAL = (
    FIXTURE_DIR
    / "random_panel_GSE239948_prestandardization.npy"
)

OUTPUT_GENES = (
    FIXTURE_DIR
    / "random_panel_background_genes.csv"
)

OUTPUT_LEGACY_BINS = (
    FIXTURE_DIR
    / "expected_legacy_random_panel_bins.csv"
)

OUTPUT_LEGACY_RANDOM = (
    FIXTURE_DIR
    / "expected_legacy_random_panel_controls.csv"
)

OUTPUT_LEGACY_CLASSIFICATION = (
    FIXTURE_DIR
    / "expected_legacy_representation_classification.csv"
)

OUTPUT_MANIFEST = (
    FIXTURE_DIR
    / "random_panel_fixture_manifest.json"
)

PRIMARY_MODULES = [
    "M34",
    "M11",
    "M24",
    "M40",
]

EXPECTED_SHARED_GENES = 13484
EXPECTED_CANDIDATE_GENES = 13212
EXPECTED_RANDOM_PANELS = 1000

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


def read_json(
    path: Path,
) -> dict[str, Any]:
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


def legacy_output_hash(
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

    for raw_path, metadata in (
        outputs.items()
    ):
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
            f"legacy hash for {basename}: "
            f"{matches}"
        )

    return matches[0]


def verify_legacy_output(
    manifest: dict[str, Any],
    path: Path,
) -> str:
    expected = legacy_output_hash(
        manifest,
        path.name,
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
    return {
        "sha256":
            sha256_file(path),
        "size_bytes":
            path.stat().st_size,
    }


def main() -> None:
    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- build authoritative random-panel fixture"
    )
    print("=" * 80)

    print(
        f"Tool version: "
        f"{TOOL_VERSION}"
    )

    for path in [
        LEGACY_RANDOM_FILE,
        LEGACY_CLASSIFICATION_FILE,
        LEGACY_MANIFEST_FILE,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Required file missing:\n"
                f"{path}"
            )

    legacy_manifest = read_json(
        LEGACY_MANIFEST_FILE
    )

    if (
        str(
            legacy_manifest.get(
                "script_version",
                "",
            )
        )
        != EXPECTED_SCRIPT_VERSION
    ):
        raise RuntimeError(
            "Unexpected legacy script version."
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

    if int(
        legacy_manifest.get(
            "random_panels",
            -1,
        )
    ) != EXPECTED_RANDOM_PANELS:
        raise RuntimeError(
            "Legacy manifest does not "
            "confirm 1000 random panels."
        )

    random_sha = verify_legacy_output(
        legacy_manifest,
        LEGACY_RANDOM_FILE,
    )

    classification_sha = (
        verify_legacy_output(
            legacy_manifest,
            LEGACY_CLASSIFICATION_FILE,
        )
    )

    print("")
    print(
        "Legacy output hashes: PASS"
    )

    print("")
    print(
        "Loading authoritative script 46 v2..."
    )

    legacy = load_legacy_module()

    if (
        legacy.SCRIPT_VERSION
        != EXPECTED_SCRIPT_VERSION
    ):
        raise RuntimeError(
            "Loaded legacy module has "
            "unexpected script version."
        )

    print(
        "Reconstructing exact raw-input "
        "analysis matrices..."
    )

    inputs = (
        reconstruct_authoritative_inputs(
            legacy
        )
    )

    reference_pre = (
        inputs[
            "reference"
        ].copy()
    )

    external_pre = (
        inputs[
            "external"
        ].copy()
    )

    reference_z = legacy.zscore_columns(
        reference_pre
    )

    external_z = legacy.zscore_columns(
        external_pre
    )

    common = (
        reference_z.columns
        .intersection(
            external_z.columns
        )
    )

    if (
        len(common)
        != EXPECTED_SHARED_GENES
    ):
        raise RuntimeError(
            "Unexpected authoritative "
            "shared-gene universe:\n"
            f"expected={EXPECTED_SHARED_GENES}\n"
            f"observed={len(common)}"
        )

    reference_fixture = (
        reference_pre.loc[
            :,
            common,
        ].copy()
    )

    external_fixture = (
        external_pre.loc[
            :,
            common,
        ].copy()
    )

    weights = (
        inputs[
            "weights"
        ].copy()
    )

    frozen_union = set(
        weights.loc[
            weights[
                "module_label"
            ]
            .astype(str)
            .isin(
                PRIMARY_MODULES
            ),
            "canine_gene_symbol",
        ]
        .astype(str)
    )

    frozen_union.discard("")

    candidate_genes = [
        str(gene)
        for gene in common
        if str(gene)
        not in frozen_union
    ]

    if (
        len(candidate_genes)
        != EXPECTED_CANDIDATE_GENES
    ):
        raise RuntimeError(
            "Unexpected candidate-gene "
            "universe:\n"
            f"expected={EXPECTED_CANDIDATE_GENES}\n"
            f"observed={len(candidate_genes)}"
        )

    legacy_bins = (
        legacy.variability_bins(
            reference_z.loc[
                :,
                common,
            ],
            external_z.loc[
                :,
                common,
            ],
        )
    )

    if (
        legacy_bins.shape[0]
        != EXPECTED_SHARED_GENES
    ):
        raise RuntimeError(
            "Unexpected legacy-bin count."
        )

    print("")
    print(
        "Selected authoritative mapping:"
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

    print("")
    print(
        "Verifying authoritative random-panel "
        "replay before freezing fixture..."
    )

    locked_structure = pd.read_csv(
        (
            SOURCE_ROOT
            / "results"
            / "tables"
            / "GSE239948_external_module_structure_preservation.csv"
        )
    )

    replay = (
        legacy.random_panel_controls(
            reference_z=reference_z,
            external_z=external_z,
            weights=weights,
            observed=locked_structure,
        )
    )

    locked_random = pd.read_csv(
        LEGACY_RANDOM_FILE
    )

    replay = (
        replay.set_index(
            "module_label"
        )
        .reindex(
            PRIMARY_MODULES
        )
        .reset_index()
    )

    locked_random = (
        locked_random.set_index(
            "module_label"
        )
        .reindex(
            PRIMARY_MODULES
        )
        .reset_index()
    )

    numeric_columns = [
        "n_module_genes",
        "n_random_panels",
        "observed_edge_spearman",
        "random_edge_median",
        "random_edge_q05",
        "random_edge_q95",
        "random_panel_empirical_p",
    ]

    for column in numeric_columns:
        observed_values = pd.to_numeric(
            replay[column],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        expected_values = pd.to_numeric(
            locked_random[column],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        if not np.allclose(
            observed_values,
            expected_values,
            rtol=1.0e-8,
            atol=1.0e-10,
            equal_nan=True,
        ):
            raise RuntimeError(
                "Authoritative replay no longer "
                "matches locked random-panel "
                f"output for column {column}."
            )

    print(
        "Authoritative replay: PASS"
    )

    FIXTURE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    # Binary NumPy serialization is intentional.
    # The legacy z-variance matching is sensitive
    # to decimal text round trips at ~1e-14.
    np.save(
        OUTPUT_REFERENCE,
        reference_fixture.to_numpy(
            dtype=np.float64
        ),
        allow_pickle=False,
    )

    np.save(
        OUTPUT_EXTERNAL,
        external_fixture.to_numpy(
            dtype=np.float64
        ),
        allow_pickle=False,
    )

    pd.DataFrame(
        {
            "gene_symbol":
                [
                    str(gene)
                    for gene in common
                ]
        }
    ).to_csv(
        OUTPUT_GENES,
        index=False,
    )

    pd.DataFrame(
        {
            "gene_symbol":
                [
                    str(gene)
                    for gene
                    in legacy_bins.index
                ],
            "legacy_variability_bin":
                legacy_bins.to_numpy(),
        }
    ).to_csv(
        OUTPUT_LEGACY_BINS,
        index=False,
    )

    shutil.copy2(
        LEGACY_RANDOM_FILE,
        OUTPUT_LEGACY_RANDOM,
    )

    shutil.copy2(
        LEGACY_CLASSIFICATION_FILE,
        OUTPUT_LEGACY_CLASSIFICATION,
    )

    if (
        sha256_file(
            OUTPUT_LEGACY_RANDOM
        )
        != random_sha
    ):
        raise RuntimeError(
            "Copied random-panel oracle "
            "hash mismatch."
        )

    if (
        sha256_file(
            OUTPUT_LEGACY_CLASSIFICATION
        )
        != classification_sha
    ):
        raise RuntimeError(
            "Copied classification oracle "
            "hash mismatch."
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
        "outcome_loaded":
            False,
        "authoritative_input_path": {
            "source":
                "raw_GSE239948_reconstruction",
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
        "contract": {
            "modules":
                PRIMARY_MODULES,
            "shared_genes":
                EXPECTED_SHARED_GENES,
            "candidate_genes":
                EXPECTED_CANDIDATE_GENES,
            "random_panels":
                EXPECTED_RANDOM_PANELS,
            "variability_bins":
                10,
            "random_seed":
                42,
            "legacy_matching_mode":
                "legacy_z_variance",
            "corrected_matching_mode":
                (
                    "prestandardization_"
                    "variance"
                ),
            "legacy_numerical_warning":
                (
                    "Legacy variability bins are "
                    "defined from variance after "
                    "gene-wise z-standardization "
                    "and are sensitive to decimal "
                    "serialization round trips."
                ),
            "fixture_serialization":
                (
                    "float64_npy_before_"
                    "z_standardization"
                ),
        },
        "dimensions": {
            "reference_samples":
                int(
                    reference_fixture.shape[0]
                ),
            "external_samples":
                int(
                    external_fixture.shape[0]
                ),
            "shared_genes":
                int(
                    len(common)
                ),
        },
        "fixture_files": {},
    }

    for path in [
        OUTPUT_REFERENCE,
        OUTPUT_EXTERNAL,
        OUTPUT_GENES,
        OUTPUT_LEGACY_BINS,
        OUTPUT_LEGACY_RANDOM,
        OUTPUT_LEGACY_CLASSIFICATION,
    ]:
        manifest[
            "fixture_files"
        ][path.name] = (
            file_record(path)
        )

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
        "Authoritative random-panel "
        "fixture summary"
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
        "Shared genes: "
        f"{len(common)}"
    )

    print(
        "Candidate genes after frozen "
        "exclusion: "
        f"{len(candidate_genes)}"
    )

    print(
        "Legacy bins frozen: "
        f"{legacy_bins.shape[0]}"
    )

    print(
        "Authoritative replay: PASS"
    )

    print(
        "Fixture serialization: "
        "float64 NPY"
    )

    print(
        "Manifest: "
        + str(
            OUTPUT_MANIFEST
            .relative_to(ROOT)
        )
    )

    print("=" * 80)


if __name__ == "__main__":
    main()
