"""File adapter for frozen-program reliability analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

import transport_audit
from transport_audit.preservation_reliability import (
    analyze_reliability,
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
        "--reliability",
        required=True,
    )
    parser.add_argument(
        "--loo",
        required=True,
    )
    parser.add_argument(
        "--manifest",
        required=True,
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    external_path = Path(
        args.external
    )

    weights_path = Path(
        args.weights
    )

    config_path = Path(
        args.config
    )

    reliability_path = Path(
        args.reliability
    )

    loo_path = Path(
        args.loo
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
            "Reliability configuration "
            "must explicitly set "
            "outcome_loaded=false."
        )

    external = pd.read_csv(
        external_path,
        index_col=0,
    )

    weights = pd.read_csv(
        weights_path
    )

    external.index = (
        external.index
        .astype(str)
    )

    modules = [
        str(module)
        for module
        in config[
            "module_order"
        ]
    ]

    reliability, loo = (
        analyze_reliability(
            external_expression=external,
            weights=weights,
            module_order=modules,
            minimum_genes=int(
                config[
                    "minimum_genes"
                ]
            ),
            n_split_half_repeats=int(
                config[
                    "split_half_repeats"
                ]
            ),
            base_seed=int(
                config[
                    "base_seed"
                ]
            ),
            module_seed_stride=int(
                config[
                    "module_seed_stride"
                ]
            ),
            split_half_seed_offset=int(
                config[
                    "split_half_seed_offset"
                ]
            ),
        )
    )

    reliability.to_csv(
        reliability_path,
        index=False,
    )

    loo.to_csv(
        loo_path,
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
                "frozen_program_"
                "reliability"
            ),
        "outcome_loaded":
            False,
        "configuration":
            config,
        "inputs": {
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
            "reliability": {
                "sha256":
                    sha256_file(
                        reliability_path
                    ),
                "rows":
                    int(
                        reliability.shape[0]
                    ),
            },
            "gene_leave_one_out": {
                "sha256":
                    sha256_file(
                        loo_path
                    ),
                "rows":
                    int(
                        loo.shape[0]
                    ),
            },
        },
        "guardrails": [
            (
                "No outcome data were "
                "loaded."
            ),
            (
                "Frozen program membership "
                "and loadings were unchanged."
            ),
            (
                "Split-half seeds were "
                "provided by the frozen "
                "configuration."
            ),
            (
                "Gene leave-one-out compares "
                "each reduced score with the "
                "full frozen signed score."
            ),
            (
                "Scientific calculations are "
                "delegated to "
                "transport_audit."
                "preservation_reliability."
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
        "Frozen preservation "
        "reliability adapter"
    )
    print("=" * 80)

    print(
        "External cohort: "
        f"{config['external_cohort']}"
    )

    print(
        "External samples: "
        f"{external.shape[0]}"
    )

    print(
        "Modules: "
        f"{reliability.shape[0]}"
    )

    print(
        "Split-half repeats: "
        f"{config['split_half_repeats']}"
    )

    print(
        "LOO rows: "
        f"{loo.shape[0]}"
    )

    print(
        "Outcome loaded: False"
    )

    print("=" * 80)


if __name__ == "__main__":
    main()
