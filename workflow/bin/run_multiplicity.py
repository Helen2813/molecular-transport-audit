from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import transport_audit
from transport_audit.multiplicity import (
    benjamini_hochberg,
    significance_mask,
)


RUNNER_VERSION = "0.1.0"


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
    parser = argparse.ArgumentParser(
        description=(
            "Apply Benjamini-Hochberg correction "
            "through the Molecular Transport Audit core."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--manifest",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--p-column",
        default="primary_p",
    )

    parser.add_argument(
        "--q-column",
        default="projectwide_q_12",
    )

    parser.add_argument(
        "--support-column",
        default="projectwide_fdr_supported",
    )

    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not args.input.exists():
        raise FileNotFoundError(
            f"Input file not found: {args.input}"
        )

    table = pd.read_csv(
        args.input
    )

    if args.p_column not in table.columns:
        raise KeyError(
            f"Missing p-value column: {args.p_column}"
        )

    pvalues = table[
        args.p_column
    ].to_numpy(
        dtype=float
    )

    qvalues = benjamini_hochberg(
        pvalues
    )

    supported = significance_mask(
        qvalues,
        alpha=args.alpha,
    )

    result = table.copy()

    result[
        args.q_column
    ] = qvalues

    result[
        args.support_column
    ] = supported

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    result.to_csv(
        args.output,
        index=False,
    )

    manifest = {
        "runner_version": RUNNER_VERSION,
        "transport_audit_version": (
            transport_audit.__version__
        ),
        "generated_at_utc": (
            datetime.now(
                timezone.utc
            ).isoformat()
        ),
        "input": {
            "name": args.input.name,
            "sha256": sha256_file(
                args.input
            ),
        },
        "output": {
            "name": args.output.name,
            "sha256": sha256_file(
                args.output
            ),
        },
        "parameters": {
            "p_column": args.p_column,
            "q_column": args.q_column,
            "support_column": (
                args.support_column
            ),
            "alpha": args.alpha,
        },
        "summary": {
            "n_rows": int(
                table.shape[0]
            ),
            "n_finite_pvalues": int(
                np.isfinite(
                    pvalues
                ).sum()
            ),
            "n_fdr_supported": int(
                supported.sum()
            ),
        },
    }

    args.manifest.write_text(
        json.dumps(
            manifest,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        f"Input rows: {table.shape[0]}"
    )

    print(
        "Finite p-values: "
        f"{np.isfinite(pvalues).sum()}"
    )

    print(
        "FDR-supported results: "
        f"{supported.sum()}"
    )

    print(
        f"Output: {args.output}"
    )


if __name__ == "__main__":
    main()
