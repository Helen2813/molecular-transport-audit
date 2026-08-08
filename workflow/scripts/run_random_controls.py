"""Nextflow adapter for corrected random-panel controls."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import transport_audit
from transport_audit.random_controls import (
    PRESTANDARDIZATION_VARIANCE,
    random_panel_controls,
)


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
        "--reference-array",
        required=True,
    )
    parser.add_argument(
        "--external-array",
        required=True,
    )
    parser.add_argument(
        "--genes",
        required=True,
    )
    parser.add_argument(
        "--weights",
        required=True,
    )
    parser.add_argument(
        "--observed-structure",
        required=True,
    )
    parser.add_argument(
        "--config",
        required=True,
    )
    parser.add_argument(
        "--output",
        required=True,
    )
    parser.add_argument(
        "--manifest",
        required=True,
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    reference_path = Path(
        args.reference_array
    )
    external_path = Path(
        args.external_array
    )
    genes_path = Path(
        args.genes
    )
    weights_path = Path(
        args.weights
    )
    structure_path = Path(
        args.observed_structure
    )
    config_path = Path(
        args.config
    )
    output_path = Path(
        args.output
    )
    manifest_path = Path(
        args.manifest
    )

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
            "Random-control configuration "
            "must set outcome_loaded=false."
        )

    matching_mode = str(
        config[
            "random_panel_matching_mode"
        ]
    )

    if (
        matching_mode
        != PRESTANDARDIZATION_VARIANCE
    ):
        raise RuntimeError(
            "Unified v0.1 audit requires "
            "prestandardization_variance."
        )

    genes = (
        pd.read_csv(
            genes_path
        )[
            "gene_symbol"
        ]
        .astype(str)
        .tolist()
    )

    reference_array = np.load(
        reference_path,
        allow_pickle=False,
    )

    external_array = np.load(
        external_path,
        allow_pickle=False,
    )

    if (
        reference_array.ndim != 2
        or external_array.ndim != 2
    ):
        raise ValueError(
            "Expression arrays must be 2D."
        )

    if (
        reference_array.shape[1]
        != len(genes)
        or external_array.shape[1]
        != len(genes)
    ):
        raise ValueError(
            "Gene list does not match "
            "expression-array dimensions."
        )

    reference = pd.DataFrame(
        reference_array,
        columns=genes,
    )

    external = pd.DataFrame(
        external_array,
        columns=genes,
    )

    weights = pd.read_csv(
        weights_path
    )

    structure = pd.read_csv(
        structure_path
    )

    result = random_panel_controls(
        reference_expression=reference,
        external_expression=external,
        weights=weights,
        observed_structure=structure,
        module_order=[
            str(module)
            for module
            in config[
                "module_order"
            ]
        ],
        matching_mode=matching_mode,
        n_random_panels=int(
            config[
                "random_panels"
            ]
        ),
        n_variability_bins=int(
            config[
                "random_panel_variability_bins"
            ]
        ),
        random_seed=int(
            config[
                "random_panel_seed"
            ]
        ),
        minimum_genes=int(
            config[
                "minimum_genes"
            ]
        ),
    )

    result.to_csv(
        output_path,
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
            "corrected_random_panel_specificity",
        "matching_mode":
            matching_mode,
        "outcome_loaded":
            False,
        "random_panels":
            int(
                config[
                    "random_panels"
                ]
            ),
        "shared_genes":
            int(
                len(genes)
            ),
        "inputs": {
            "reference_array_sha256":
                sha256_file(
                    reference_path
                ),
            "external_array_sha256":
                sha256_file(
                    external_path
                ),
            "genes_sha256":
                sha256_file(
                    genes_path
                ),
            "weights_sha256":
                sha256_file(
                    weights_path
                ),
            "observed_structure_sha256":
                sha256_file(
                    structure_path
                ),
        },
        "output": {
            "sha256":
                sha256_file(
                    output_path
                ),
            "rows":
                int(
                    result.shape[0]
                ),
        },
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
        "Corrected random-panel controls"
    )
    print("=" * 80)
    print(
        f"Shared genes: {len(genes)}"
    )
    print(
        "Random panels per module: "
        f"{config['random_panels']}"
    )
    print(
        f"Matching mode: {matching_mode}"
    )
    print(
        "Outcome loaded: False"
    )
    print("=" * 80)


if __name__ == "__main__":
    main()
