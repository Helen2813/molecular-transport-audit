from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_VERSION = "05c1-audit-classical-bootstrap-degeneracy-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

V2_DIR = DATA_ROOT / "paper4_tcbb_primary_bootstrap_reliability_v2"
CHECKPOINT = (
    V2_DIR
    / "per_program"
    / "SCANB_GSE96058"
    / "TCGA_M001_bootstrap_checkpoint_v2.npz"
)
DEGENERATE_COUNTS = (
    V2_DIR
    / "per_program"
    / "SCANB_GSE96058"
    / "TCGA_M001_bootstrap_degenerate_gene_counts_v2.tsv"
)

GENE_FILE = (
    DATA_ROOT
    / "paper4_tcbb_primary_pooled_direct_preservation_v1"
    / "per_program"
    / "SCANB_GSE96058"
    / "TCGA_M001_evaluable_genes_v1.tsv"
)
SCANB_EXPR = (
    DATA_ROOT
    / "SCANB_GSE96058"
    / "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
)
SCANB_PRIMARY = (
    DATA_ROOT
    / "paper4_tcbb_scanb_pairing_contract_v1"
    / "scanb_primary_profiles_frozen_v1.tsv"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_classical_bootstrap_degeneracy_audit_v1"


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def canon_symbol(x: object) -> str:
    if x is None:
        return ""
    return " ".join(str(x).strip().split()).upper()


def logsumexp(vals: list[float]) -> float:
    if not vals:
        return float("-inf")
    m = max(vals)
    return m + math.log(sum(math.exp(v - m) for v in vals))


def exact_constant_bootstrap_probability(values: np.ndarray) -> float:
    """
    For an ordinary size-n nonparametric bootstrap of a single gene, the
    resample is constant iff every draw comes from one exact observed value.
    P(constant) = sum_j p_j^n over distinct empirical values.
    """
    vals = np.asarray(values)
    n = len(vals)
    _, counts = np.unique(vals, return_counts=True)
    logs = [
        n * math.log(int(c) / n)
        for c in counts
        if c > 0
    ]
    return math.exp(logsumexp(logs))


def main() -> None:
    print("=" * 136)
    print("Paper 4 / TCBB - audit why the fixed-gene ordinary target bootstrap is degenerate")
    print("=" * 136)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target outcomes/treatment loaded:                    NO")
    print("  New preservation statistic calculated:               NO")
    print("  Operation: explain the 05c v2 bootstrap failure only")
    print("=" * 136)

    for p in [
        CHECKPOINT,
        DEGENERATE_COUNTS,
        GENE_FILE,
        SCANB_EXPR,
        SCANB_PRIMARY,
    ]:
        require(p)

    ck = np.load(CHECKPOINT)
    completed = int(ck["completed"])
    attempts = int(ck["attempts"])
    invalid = int(ck["invalid_attempts"])

    if completed != 0 or attempts != 10_000 or invalid != 10_000:
        raise RuntimeError(
            f"Unexpected 05c v2 failure state: completed={completed}, "
            f"attempts={attempts}, invalid={invalid}"
        )

    print(f"\n05c v2 failure replay: valid={completed}, attempts={attempts}, invalid={invalid}")

    deg = pd.read_csv(DEGENERATE_COUNTS, sep="\t")
    genes = pd.read_csv(GENE_FILE, sep="\t", dtype=str).fillna("")
    genes["Hugo_Symbol"] = genes["Hugo_Symbol"].map(canon_symbol)

    if len(genes) != 1697:
        raise RuntimeError(f"Expected 1,697 M001 evaluable genes, found {len(genes)}.")

    primary = pd.read_csv(SCANB_PRIMARY, sep="\t", dtype=str).fillna("")
    primary_titles = list(primary["primary_title"])
    if len(primary_titles) != 3273:
        raise RuntimeError(f"Expected 3,273 SCAN-B primaries, found {len(primary_titles)}.")

    wanted = set(genes["Hugo_Symbol"])
    expression_rows = []

    for chunk in pd.read_csv(
        SCANB_EXPR,
        compression="gzip",
        chunksize=1500,
        low_memory=False,
    ):
        symbol_col = chunk.columns[0]
        symbols = chunk[symbol_col].map(canon_symbol)
        mask = symbols.isin(wanted)
        if not mask.any():
            continue

        sub = chunk.loc[mask, [symbol_col] + primary_titles].copy()
        arr = sub[primary_titles].to_numpy(dtype=np.float64)

        fpkm = np.maximum(np.exp2(arr) - 0.1, 0.0)
        transformed = np.log2(fpkm + 1.0)

        tmp = pd.DataFrame(transformed, columns=primary_titles)
        tmp.insert(0, "Hugo_Symbol", symbols.loc[mask].to_numpy())
        expression_rows.append(tmp)

    if not expression_rows:
        raise RuntimeError("No M001 genes found in SCAN-B expression.")

    expr = pd.concat(expression_rows, ignore_index=True)
    expr = expr.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    missing = sorted(wanted - set(expr.index))
    if missing:
        raise RuntimeError(f"M001 expression genes missing: {missing[:20]}")

    rows = []
    for gene in genes["Hugo_Symbol"]:
        v = expr.loc[gene].to_numpy(dtype=np.float64)
        unique, counts = np.unique(v, return_counts=True)
        max_count = int(counts.max())
        max_fraction = max_count / len(v)
        p_const = exact_constant_bootstrap_probability(v)

        rows.append(
            {
                "Hugo_Symbol": gene,
                "n_samples": len(v),
                "n_unique_exact_values": int(len(unique)),
                "most_common_exact_value_count": max_count,
                "most_common_exact_value_fraction": max_fraction,
                "ordinary_bootstrap_constant_probability": p_const,
            }
        )

    diag = pd.DataFrame(rows)
    deg["Hugo_Symbol"] = deg["Hugo_Symbol"].map(canon_symbol)
    diag = diag.merge(
        deg[["Hugo_Symbol", "degenerate_attempt_count"]],
        on="Hugo_Symbol",
        how="left",
        validate="one_to_one",
    )
    diag["degenerate_attempt_count"] = (
        diag["degenerate_attempt_count"].fillna(0).astype(int)
    )
    diag["observed_degenerate_attempt_fraction"] = (
        diag["degenerate_attempt_count"] / attempts
    )

    diag = diag.sort_values(
        ["degenerate_attempt_count", "ordinary_bootstrap_constant_probability"],
        ascending=False,
    ).reset_index(drop=True)

    expected_degenerate_genes_per_draw = float(
        diag["ordinary_bootstrap_constant_probability"].sum()
    )
    genes_ever = int((diag["degenerate_attempt_count"] > 0).sum())
    genes_always = int((diag["degenerate_attempt_count"] == attempts).sum())

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    detail_out = OUT_DIR / "scanb_m001_classical_bootstrap_degeneracy_gene_audit_v1.tsv"
    diag.to_csv(detail_out, sep="\t", index=False)

    summary = {
        "script_version": SCRIPT_VERSION,
        "status": "CLASSICAL_FIXED_GENE_BOOTSTRAP_DEGENERACY_CONFIRMED",
        "target": "SCANB_GSE96058",
        "program_id": "TCGA_M001",
        "n_samples": 3273,
        "n_fixed_genes": 1697,
        "classical_bootstrap_attempts": attempts,
        "valid_attempts": completed,
        "invalid_attempts": invalid,
        "genes_ever_degenerate": genes_ever,
        "genes_degenerate_in_all_attempts": genes_always,
        "sum_of_per_gene_exact_constant_probabilities": expected_degenerate_genes_per_draw,
        "interpretation": (
            "The fixed-gene ordinary nonparametric bootstrap is not practically "
            "estimable for this high-dimensional sparse-expression program. "
            "No preservation effect or interval was used to make this determination."
        ),
    }

    summary_out = OUT_DIR / "classical_bootstrap_degeneracy_audit_v1.json"
    summary_out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print()
    print("Top genes by observed degenerate-attempt count:")
    for r in diag.head(15).itertuples(index=False):
        print(
            f"  {r.Hugo_Symbol}: "
            f"degenerate={r.degenerate_attempt_count}/{attempts}; "
            f"max-value fraction={r.most_common_exact_value_fraction:.6f}; "
            f"single-gene P(constant bootstrap)={r.ordinary_bootstrap_constant_probability:.6g}"
        )

    print()
    print("=" * 136)
    print("05c1 CLASSICAL BOOTSTRAP DEGENERACY AUDIT: COMPLETE")
    print("=" * 136)
    print(f"Genes ever degenerate:          {genes_ever}/{len(diag)}")
    print(f"Genes degenerate in all draws: {genes_always}/{len(diag)}")
    print(
        "Sum of per-gene exact P(constant): "
        f"{expected_degenerate_genes_per_draw:.4f}"
    )
    print("No preservation effect or uncertainty interval was calculated.")
    print()
    print(f"Detail:  {detail_out}")
    print(f"Summary: {summary_out}")
    print("=" * 136)


if __name__ == "__main__":
    main()
