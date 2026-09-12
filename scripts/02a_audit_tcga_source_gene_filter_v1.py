from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_VERSION = "02a-audit-tcga-source-gene-filter-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")
TCGA_EXPR = (
    DATA_ROOT
    / "paper4_tcbb_input_audit_v1"
    / "staged_continuous_inputs"
    / "TCGA_BRCA_PanCanAtlas2018"
    / "data_mrna_seq_v2_rsem.txt"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_tcga_source_filter_audit_v1"

CHUNK_ROWS = 256


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def robust_mad(x: np.ndarray) -> np.ndarray:
    """
    Raw median absolute deviation per row (no 1.4826 scale factor).
    x: genes x samples
    """
    med = np.nanmedian(x, axis=1)
    return np.nanmedian(np.abs(x - med[:, None]), axis=1)


def qdict(x: np.ndarray) -> dict[str, float | None]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {
            "min": None, "q01": None, "q05": None, "q25": None,
            "median": None, "q75": None, "q95": None, "q99": None, "max": None
        }
    q = np.quantile(x, [0, .01, .05, .25, .50, .75, .95, .99, 1])
    keys = ["min", "q01", "q05", "q25", "median", "q75", "q95", "q99", "max"]
    return {k: float(v) for k, v in zip(keys, q)}


def banner() -> None:
    print("=" * 116)
    print("Paper 4 / TCBB - audit TCGA-BRCA source-only gene-filter candidates")
    print("=" * 116)
    print(f"Script version: {SCRIPT_VERSION}")
    print(f"TCGA source:    {TCGA_EXPR}")
    print()
    print("Scientific guard:")
    print("  SCAN-B opened:                                      NO")
    print("  METABRIC opened:                                    NO")
    print("  GSE239948 opened:                                   NO")
    print("  Clinical outcomes loaded:                           NO")
    print("  Target preservation statistics calculated:          NO")
    print("  Target coverage used to select genes:                NO")
    print("  Operation: source-only expression/measurability audit")
    print("=" * 116)


def main() -> None:
    banner()
    require(TCGA_EXPR)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    header = pd.read_csv(TCGA_EXPR, sep="\t", nrows=0, low_memory=False)
    cols = list(header.columns)
    if cols[:2] != ["Hugo_Symbol", "Entrez_Gene_Id"]:
        raise RuntimeError(f"Unexpected TCGA identifier columns: {cols[:2]}")
    sample_cols = cols[2:]
    n_samples = len(sample_cols)

    rows = []

    for chunk_i, df in enumerate(
        pd.read_csv(TCGA_EXPR, sep="\t", chunksize=CHUNK_ROWS, low_memory=False),
        start=1,
    ):
        symbols = df["Hugo_Symbol"].astype("string").fillna("").str.strip()
        entrez = df["Entrez_Gene_Id"].astype("string").fillna("").str.strip()

        raw = df[sample_cols].apply(pd.to_numeric, errors="coerce").to_numpy(np.float64)
        finite = np.isfinite(raw)

        if np.nanmin(raw) < 0:
            raise RuntimeError("Negative TCGA RSEM value encountered; expression contract violated.")

        logx = np.log2(raw + 1.0)

        finite_n = finite.sum(axis=1)
        missing_fraction = 1.0 - (finite_n / n_samples)

        # Source-only measurability summaries on raw RSEM scale.
        frac_gt0 = np.nansum(raw > 0.0, axis=1) / n_samples
        frac_ge1 = np.nansum(raw >= 1.0, axis=1) / n_samples
        frac_ge5 = np.nansum(raw >= 5.0, axis=1) / n_samples
        frac_ge10 = np.nansum(raw >= 10.0, axis=1) / n_samples

        with np.errstate(invalid="ignore", divide="ignore"):
            mean_log = np.nanmean(logx, axis=1)
            median_log = np.nanmedian(logx, axis=1)
            var_log = np.nanvar(logx, axis=1, ddof=1)
            q25 = np.nanpercentile(logx, 25, axis=1)
            q75 = np.nanpercentile(logx, 75, axis=1)
            iqr_log = q75 - q25
            mad_log = robust_mad(logx)

        for i in range(len(df)):
            rows.append(
                {
                    "Hugo_Symbol": str(symbols.iloc[i]),
                    "Entrez_Gene_Id": str(entrez.iloc[i]),
                    "missing_fraction": float(missing_fraction[i]),
                    "frac_rsem_gt0": float(frac_gt0[i]),
                    "frac_rsem_ge1": float(frac_ge1[i]),
                    "frac_rsem_ge5": float(frac_ge5[i]),
                    "frac_rsem_ge10": float(frac_ge10[i]),
                    "mean_log2_rsem1": float(mean_log[i]),
                    "median_log2_rsem1": float(median_log[i]),
                    "mad_log2_rsem1": float(mad_log[i]),
                    "iqr_log2_rsem1": float(iqr_log[i]),
                    "var_log2_rsem1": float(var_log[i]),
                }
            )

        if chunk_i % 25 == 0:
            print(f"  audited {min(chunk_i * CHUNK_ROWS, 20531):,} genes ...")

    metrics = pd.DataFrame(rows)

    # Identifier integrity.
    blank_symbol = metrics["Hugo_Symbol"].eq("")
    duplicate_symbol = metrics["Hugo_Symbol"].duplicated(keep=False) & ~blank_symbol
    duplicate_entrez = (
        metrics["Entrez_Gene_Id"].ne("")
        & metrics["Entrez_Gene_Id"].duplicated(keep=False)
    )

    # Candidate source-only filters. These are AUDITED, not yet frozen.
    finite_ok = metrics["missing_fraction"] <= 0.05
    nonzero_mad = metrics["mad_log2_rsem1"] > 0

    candidate_defs = {
        "C0_finite_and_nonzero_MAD": finite_ok & nonzero_mad,
        "C1_RSEM_ge1_in_10pct_and_nonzero_MAD": (
            finite_ok & (metrics["frac_rsem_ge1"] >= 0.10) & nonzero_mad
        ),
        "C2_RSEM_ge1_in_20pct_and_nonzero_MAD": (
            finite_ok & (metrics["frac_rsem_ge1"] >= 0.20) & nonzero_mad
        ),
        "C3_RSEM_ge1_in_50pct_and_nonzero_MAD": (
            finite_ok & (metrics["frac_rsem_ge1"] >= 0.50) & nonzero_mad
        ),
        "C4_RSEM_ge5_in_20pct_and_nonzero_MAD": (
            finite_ok & (metrics["frac_rsem_ge5"] >= 0.20) & nonzero_mad
        ),
        "C5_RSEM_ge10_in_20pct_and_nonzero_MAD": (
            finite_ok & (metrics["frac_rsem_ge10"] >= 0.20) & nonzero_mad
        ),
    }

    # Add robust-variance ranking variants only within C2, so we can see whether
    # a fixed top-MAD cap is even needed computationally.
    c2_mask = candidate_defs["C2_RSEM_ge1_in_20pct_and_nonzero_MAD"]
    c2 = metrics.loc[c2_mask].copy()
    c2_ranked = c2.sort_values(
        ["mad_log2_rsem1", "var_log2_rsem1", "Hugo_Symbol"],
        ascending=[False, False, True],
        kind="mergesort",
    )

    summary_rows = []
    for name, mask in candidate_defs.items():
        subset = metrics.loc[mask]
        summary_rows.append(
            {
                "candidate": name,
                "n_genes": int(mask.sum()),
                "fraction_of_20531": float(mask.mean()),
                "median_frac_rsem_ge1": (
                    float(subset["frac_rsem_ge1"].median()) if len(subset) else None
                ),
                "median_MAD_log2": (
                    float(subset["mad_log2_rsem1"].median()) if len(subset) else None
                ),
                "median_variance_log2": (
                    float(subset["var_log2_rsem1"].median()) if len(subset) else None
                ),
            }
        )

    for cap in [5000, 7500, 10000, 12500, 15000]:
        n = min(cap, len(c2_ranked))
        if n == 0:
            continue
        subset = c2_ranked.iloc[:n]
        summary_rows.append(
            {
                "candidate": f"C2_plus_top_{cap}_by_MAD",
                "n_genes": int(n),
                "fraction_of_20531": float(n / len(metrics)),
                "median_frac_rsem_ge1": float(subset["frac_rsem_ge1"].median()),
                "median_MAD_log2": float(subset["mad_log2_rsem1"].median()),
                "median_variance_log2": float(subset["var_log2_rsem1"].median()),
            }
        )

    summary = pd.DataFrame(summary_rows)

    metrics_out = OUT_DIR / "tcga_source_gene_metrics_v1.tsv"
    summary_out = OUT_DIR / "tcga_source_filter_candidate_summary_v1.tsv"
    c2_ranked_out = OUT_DIR / "tcga_C2_genes_ranked_by_MAD_v1.tsv"

    metrics.to_csv(metrics_out, sep="\t", index=False)
    summary.to_csv(summary_out, sep="\t", index=False)
    c2_ranked.to_csv(c2_ranked_out, sep="\t", index=False)

    audit = {
        "script_version": SCRIPT_VERSION,
        "status": "SOURCE_ONLY_FILTER_AUDIT_NOT_YET_FROZEN",
        "scientific_guard": {
            "targets_opened": False,
            "outcomes_loaded": False,
            "target_coverage_used_for_selection": False,
            "preservation_statistics_calculated": False,
        },
        "tcga": {
            "genes": int(len(metrics)),
            "samples": int(n_samples),
            "blank_Hugo_symbols": int(blank_symbol.sum()),
            "rows_in_duplicate_Hugo_symbol_groups": int(duplicate_symbol.sum()),
            "rows_in_duplicate_Entrez_groups": int(duplicate_entrez.sum()),
            "missing_fraction_quantiles": qdict(metrics["missing_fraction"].to_numpy()),
            "frac_rsem_ge1_quantiles": qdict(metrics["frac_rsem_ge1"].to_numpy()),
            "mad_log2_quantiles": qdict(metrics["mad_log2_rsem1"].to_numpy()),
            "variance_log2_quantiles": qdict(metrics["var_log2_rsem1"].to_numpy()),
        },
        "candidate_filter_counts": {
            row["candidate"]: int(row["n_genes"])
            for row in summary_rows
        },
        "recommended_next_decision": (
            "Choose and freeze ONE source-only measurability filter before "
            "WGCNA soft-threshold diagnostics. Do not use target coverage or "
            "target preservation to choose the filter."
        ),
    }

    json_out = OUT_DIR / "tcga_source_filter_audit_v1.json"
    json_out.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print("\n" + "=" * 116)
    print("02a TCGA SOURCE-ONLY FILTER AUDIT: PASS")
    print("=" * 116)
    print(f"TCGA genes:   {len(metrics):,}")
    print(f"TCGA samples: {n_samples:,}")
    print(f"Blank Hugo symbols:                     {int(blank_symbol.sum()):,}")
    print(f"Rows in duplicate Hugo-symbol groups:   {int(duplicate_symbol.sum()):,}")
    print(f"Rows in duplicate Entrez-ID groups:     {int(duplicate_entrez.sum()):,}")
    print()
    print("Candidate source-only filters:")
    for _, row in summary.iterrows():
        print(
            f"  {row['candidate']:<43} "
            f"n={int(row['n_genes']):>6,}  "
            f"median_MAD={row['median_MAD_log2']:.6g}"
        )
    print()
    print("No filter has been frozen yet.")
    print("No target dataset was opened.")
    print("No WGCNA modules or preservation statistics were calculated.")
    print()
    print("Outputs:")
    print(f"  {summary_out}")
    print(f"  {metrics_out}")
    print(f"  {c2_ranked_out}")
    print(f"  {json_out}")
    print("=" * 116)


if __name__ == "__main__":
    main()
