"""Nextflow adapter for preservation permutation inference."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

import transport_audit
from transport_audit.preservation_inference import (
    analyze_permutation_preservation,
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

    reference_path = Path(args.reference)
    external_path = Path(args.external)
    weights_path = Path(args.weights)
    config_path = Path(args.config)
    output_path = Path(args.output)
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
            "Permutation configuration "
            "must set outcome_loaded=false."
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
        in config[
            "module_order"
        ]
    ]

    result = (
        analyze_permutation_preservation(
            reference_expression=reference,
            external_expression=external,
            weights=weights,
            module_order=modules,
            minimum_genes=int(
                config[
                    "minimum_genes"
                ]
            ),
            n_permutations=int(
                config[
                    "permutations"
                ]
            ),
            base_seed=int(
                config[
                    "permutation_base_seed"
                ]
            ),
            module_seed_stride=int(
                config[
                    "module_seed_stride"
                ]
            ),
            pca_random_state=int(
                config[
                    "pca_random_state"
                ]
            ),
        )
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
            "preservation_permutation_inference",
        "outcome_loaded":
            False,
        "modules":
            modules,
        "n_permutations":
            int(
                config[
                    "permutations"
                ]
            ),
        "inputs": {
            "reference_sha256":
                sha256_file(
                    reference_path
                ),
            "external_sha256":
                sha256_file(
                    external_path
                ),
            "weights_sha256":
                sha256_file(
                    weights_path
                ),
            "config_sha256":
                sha256_file(
                    config_path
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
        "Preservation permutation inference"
    )
    print("=" * 80)
    print(
        f"Modules: {len(modules)}"
    )
    print(
        "Permutations per test: "
        f"{config['permutations']}"
    )
    print(
        "Outcome loaded: False"
    )
    print("=" * 80)


if __name__ == "__main__":
    main()
