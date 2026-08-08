from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

import transport_audit
from transport_audit.scoring import score_programs
from transport_audit.schemas import ScoringRule


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
        "--expression",
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
        "--scores",
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

    expression_path = Path(
        args.expression
    )
    weights_path = Path(
        args.weights
    )
    config_path = Path(
        args.config
    )

    scores_path = Path(
        args.scores
    )
    coverage_path = Path(
        args.coverage
    )
    manifest_path = Path(
        args.manifest
    )

    config = json.loads(
        config_path.read_text(
            encoding="utf-8"
        )
    )

    expression = pd.read_csv(
        expression_path,
        index_col=0,
    )

    weights = pd.read_csv(
        weights_path
    )

    expression.index = (
        expression.index
        .astype(str)
    )

    module_order = [
        str(module)
        for module
        in config[
            "module_order"
        ]
    ]

    rule = ScoringRule(
        mapping_name=str(
            config[
                "mapping_name"
            ]
        ),
        minimum_genes=int(
            config[
                "minimum_genes"
            ]
        ),
        minimum_fraction=float(
            config[
                "minimum_fraction"
            ]
        ),
    )

    scores, coverage = score_programs(
        expression=expression,
        weights=weights,
        cohort=str(
            config["cohort"]
        ),
        module_order=module_order,
        rule=rule,
    )

    scores.to_csv(
        scores_path,
        index=True,
        index_label=(
            scores.index.name
            or "sample_id"
        ),
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
            "frozen_program_scoring",
        "outcome_loaded":
            False,
        "configuration":
            config,
        "input": {
            "expression": {
                "sha256":
                    sha256_file(
                        expression_path
                    ),
                "n_samples":
                    int(
                        expression.shape[0]
                    ),
                "n_genes":
                    int(
                        expression.shape[1]
                    ),
            },
            "weights": {
                "sha256":
                    sha256_file(
                        weights_path
                    ),
                "n_rows":
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
        "output": {
            "scores": {
                "sha256":
                    sha256_file(
                        scores_path
                    ),
                "n_samples":
                    int(
                        scores.shape[0]
                    ),
                "n_score_columns":
                    int(
                        scores.shape[1]
                    ),
                "columns":
                    list(
                        scores.columns
                    ),
            },
            "coverage": {
                "sha256":
                    sha256_file(
                        coverage_path
                    ),
                "n_rows":
                    int(
                        coverage.shape[0]
                    ),
            },
        },
        "guardrails": [
            (
                "Frozen program membership "
                "and loadings are unchanged."
            ),
            (
                "No outcome data are loaded."
            ),
            (
                "Coverage decisions use only "
                "expression availability."
            ),
            (
                "The adapter delegates all "
                "scientific scoring logic to "
                "transport_audit.scoring."
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
        "Frozen program scoring adapter"
    )
    print("=" * 80)
    print(
        f"Cohort: "
        f"{config['cohort']}"
    )
    print(
        f"Samples: "
        f"{scores.shape[0]}"
    )
    print(
        f"Score columns: "
        f"{scores.shape[1]}"
    )
    print(
        f"Coverage rows: "
        f"{coverage.shape[0]}"
    )
    print(
        "Outcome loaded: False"
    )
    print("=" * 80)


if __name__ == "__main__":
    main()
