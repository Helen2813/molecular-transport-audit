"""File adapter for deterministic molecular-program preservation."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

import transport_audit
from transport_audit.preservation import analyze_preservation


SCRIPT_VERSION = "0.1.0"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--reference",
        required=True,
    )
    parser.add_argument(
        "--external",
        required=True,
    )
    parser.add_argument(
        "--weights",
        required=True,
    )
    parser.add_argument(
        "--config",
        required=True,
    )
    parser.add_argument(
        "--structure",
        required=True,
    )
    parser.add_argument(
        "--coverage",
        required=True,
    )
    parser.add_argument(
        "--manifest",
        required=True,
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    reference_path = Path(args.reference)
    external_path = Path(args.external)
    weights_path = Path(args.weights)
    config_path = Path(args.config)

    structure_path = Path(args.structure)
    coverage_path = Path(args.coverage)
    manifest_path = Path(args.manifest)

    config = json.loads(
        config_path.read_text(
            encoding="utf-8"
        )
    )

    if bool(
        config.get(
            "outcome_loaded",
            True,
        )
    ):
        raise RuntimeError(
            "Preservation configuration "
            "must explicitly set "
            "outcome_loaded=false."
        )

    reference = pd.read_csv(
        reference_path,
        index_col=0,
    )

    external = pd.read_csv(
        external_path,
        index_col=0,
    )

    weights = pd.read_csv(
        weights_path
    )

    reference.index = (
        reference.index.astype(str)
    )

    external.index = (
        external.index.astype(str)
    )

    modules = [
        str(module)
        for module
        in config["module_order"]
    ]

    structure, coverage = (
        analyze_preservation(
            reference_expression=reference,
            external_expression=external,
            weights=weights,
            module_order=modules,
            minimum_genes=int(
                config["minimum_genes"]
            ),
            pca_random_state=int(
                config[
                    "pca_random_state"
                ]
            ),
        )
    )

    structure.to_csv(
        structure_path,
        index=False,
    )

    coverage.to_csv(
        coverage_path,
        index=False,
    )

    manifest = {
        "script_version":
            SCRIPT_VERSION,
        "transport_audit_version":
            transport_audit.__version__,
        "created_at_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),
        "scientific_role":
            (
                "deterministic_direct_"
                "molecular_representation_"
                "preservation"
            ),
        "outcome_loaded":
            False,
        "configuration":
            config,
        "inputs": {
            "reference": {
                "sha256":
                    sha256_file(
                        reference_path
                    ),
                "samples":
                    int(
                        reference.shape[0]
                    ),
                "genes":
                    int(
                        reference.shape[1]
                    ),
            },
            "external": {
                "sha256":
                    sha256_file(
                        external_path
                    ),
                "samples":
                    int(
                        external.shape[0]
                    ),
                "genes":
                    int(
                        external.shape[1]
                    ),
            },
            "weights": {
                "sha256":
                    sha256_file(
                        weights_path
                    ),
                "rows":
                    int(
                        weights.shape[0]
                    ),
            },
            "config": {
                "sha256":
                    sha256_file(
                        config_path
                    ),
            },
        },
        "outputs": {
            "structure": {
                "sha256":
                    sha256_file(
                        structure_path
                    ),
                "rows":
                    int(
                        structure.shape[0]
                    ),
            },
            "coverage": {
                "sha256":
                    sha256_file(
                        coverage_path
                    ),
                "rows":
                    int(
                        coverage.shape[0]
                    ),
            },
        },
        "guardrails": [
            (
                "No outcome data were "
                "loaded."
            ),
            (
                "Frozen module membership "
                "and loadings were not "
                "modified."
            ),
            (
                "External PC1 orientation "
                "used only the frozen "
                "signed score."
            ),
            (
                "This process contains no "
                "permutation, split-half, "
                "or random-panel inference."
            ),
            (
                "Scientific calculations "
                "are delegated to "
                "transport_audit.preservation."
            ),
        ],
    }

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("=" * 80)
    print(
        "Deterministic preservation adapter"
    )
    print("=" * 80)
    print(
        "Reference cohort: "
        f"{config['reference_cohort']}"
    )
    print(
        "External cohort: "
        f"{config['external_cohort']}"
    )
    print(
        "Reference samples: "
        f"{reference.shape[0]}"
    )
    print(
        "External samples: "
        f"{external.shape[0]}"
    )
    print(
        "Modules: "
        f"{structure.shape[0]}"
    )
    print(
        "Outcome loaded: False"
    )
    print("=" * 80)


if __name__ == "__main__":
    main()
