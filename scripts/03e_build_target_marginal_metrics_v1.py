from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_VERSION = "03e-build-target-marginal-metrics-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

SOURCE_METRICS = (
    DATA_ROOT
    / "paper4_tcbb_tcga_source_filter_audit_v1"
    / "tcga_source_gene_metrics_v1.tsv"
)
SOURCE_UNIVERSE = (
    DATA_ROOT
    / "paper4_tcbb_tcga_source_universe_v1"
    / "tcga_source_gene_universe_frozen_v1.tsv"
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
SCANB_MAP = (
    DATA_ROOT
    / "paper4_tcbb_target_mapping_v1"
    / "scanb_gene_mapping_frozen_v1.tsv"
)

METABRIC_EXPR = (
    DATA_ROOT
    / "paper4_tcbb_input_audit_v1"
    / "staged_continuous_inputs"
    / "METABRIC"
    / "data_mrna_illumina_microarray.txt"
)
METABRIC_MAP = (
    DATA_ROOT
    / "paper4_tcbb_target_mapping_v1"
    / "metabric_gene_mapping_frozen_v1.tsv"
)

NULL_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_mapping_null_contract_v1"
    / "structure_preserving_mapping_null_contract_v1.json"
)

OUT_DIR = DATA_ROOT / "paper4_tcbb_target_marginal_metrics_v1"

CHUNK_ROWS = 128


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def raw_mad(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    med = np.median(x)
    return float(np.median(np.abs(x - med)))


def pct_rank_average(x: pd.Series) -> pd.Series:
    return x.rank(method="average", pct=True)


def quantile_bin(score: pd.Series, q: int = 5) -> pd.Series:
    # Rank first to make qcut deterministic in the presence of many ties.
    ranks = score.rank(method="first")
    return pd.qcut(ranks, q=q, labels=False, duplicates="raise").astype(int) + 1


def prepare_source_features() -> pd.DataFrame:
    metrics = pd.read_csv(SOURCE_METRICS, sep="\t", dtype={"Hugo_Symbol": str})
    universe = pd.read_csv(SOURCE_UNIVERSE, sep="\t", dtype={"Hugo_Symbol": str})

    wanted = universe[["Hugo_Symbol"]].copy()
    wanted["Hugo_Symbol"] = wanted["Hugo_Symbol"].astype(str).str.upper().str.strip()

    metrics["Hugo_Symbol"] = metrics["Hugo_Symbol"].fillna("").astype(str).str.upper().str.strip()

    # Source universe was frozen by source row; one symbol each.
    merged = wanted.merge(
        metrics[["Hugo_Symbol", "mean_log2_rsem1", "mad_log2_rsem1"]],
        on="Hugo_Symbol",
        how="left",
        validate="one_to_one",
    )
    if merged[["mean_log2_rsem1", "mad_log2_rsem1"]].isna().any().any():
        raise RuntimeError("Missing source marginal metric in frozen 10k universe.")
    return merged


def scanb_metrics() -> pd.DataFrame:
    primary = pd.read_csv(SCANB_PRIMARY, sep="\t", dtype=str)
    primary_cols = primary["primary_title"].tolist()

    mapping = pd.read_csv(SCANB_MAP, sep="\t", dtype=str)
    row_to_symbol = {}
    dup_symbols = set()
    for _, r in mapping.iterrows():
        sym = str(r["target_normalized_symbol"]).upper().strip()
        idxs = [int(x) for x in str(r["target_row_indices_0based"]).split(",")]
        for idx in idxs:
            row_to_symbol[idx] = sym
        if len(idxs) > 1:
            dup_symbols.add(sym)

    header = pd.read_csv(SCANB_EXPR, compression="gzip", nrows=0)
    all_cols = list(header.columns)
    id_col = all_cols[0]
    missing = [c for c in primary_cols if c not in all_cols]
    if missing:
        raise RuntimeError(f"SCAN-B primary columns missing: {missing[:10]}")

    rows_out = []
    dup_vectors = defaultdict(list)
    global_row = 0

    usecols = [id_col] + primary_cols
    for chunk in pd.read_csv(
        SCANB_EXPR,
        compression="gzip",
        usecols=usecols,
        chunksize=CHUNK_ROWS,
        low_memory=False,
    ):
        for _, row in chunk.iterrows():
            if global_row in row_to_symbol:
                sym = row_to_symbol[global_row]
                y = pd.to_numeric(row[primary_cols], errors="coerce").to_numpy(float)
                fpkm = np.maximum(np.exp2(y) - 0.1, 0.0)
                tx = np.log2(fpkm + 1.0)
                if sym in dup_symbols:
                    dup_vectors[sym].append(tx)
                else:
                    finite = tx[np.isfinite(tx)]
                    rows_out.append(
                        {
                            "target_normalized_symbol": sym,
                            "target_mean": float(np.mean(finite)),
                            "target_mad": raw_mad(finite),
                            "target_missing_fraction": float(1 - finite.size / tx.size),
                        }
                    )
            global_row += 1

    for sym, vecs in dup_vectors.items():
        mat = np.vstack(vecs)
        avg = np.nanmean(mat, axis=0)
        finite = avg[np.isfinite(avg)]
        rows_out.append(
            {
                "target_normalized_symbol": sym,
                "target_mean": float(np.mean(finite)),
                "target_mad": raw_mad(finite),
                "target_missing_fraction": float(1 - finite.size / avg.size),
            }
        )

    out = pd.DataFrame(rows_out)
    if out["target_normalized_symbol"].duplicated().any():
        raise RuntimeError("SCAN-B marginal table still has duplicate symbols.")
    return out


def metabric_metrics() -> pd.DataFrame:
    mapping = pd.read_csv(METABRIC_MAP, sep="\t", dtype=str)
    row_to_symbol = {}
    dup_symbols = set()
    for _, r in mapping.iterrows():
        sym = str(r["target_normalized_symbol"]).upper().strip()
        idxs = [int(x) for x in str(r["target_row_indices_0based"]).split(",")]
        for idx in idxs:
            row_to_symbol[idx] = sym
        if len(idxs) > 1:
            dup_symbols.add(sym)

    header = pd.read_csv(METABRIC_EXPR, sep="\t", nrows=0)
    all_cols = list(header.columns)
    sample_cols = all_cols[2:]

    rows_out = []
    dup_vectors = defaultdict(list)
    global_row = 0

    for chunk in pd.read_csv(
        METABRIC_EXPR,
        sep="\t",
        chunksize=CHUNK_ROWS,
        low_memory=False,
    ):
        for _, row in chunk.iterrows():
            if global_row in row_to_symbol:
                sym = row_to_symbol[global_row]
                tx = pd.to_numeric(row[sample_cols], errors="coerce").to_numpy(float)
                if sym in dup_symbols:
                    dup_vectors[sym].append(tx)
                else:
                    finite = tx[np.isfinite(tx)]
                    rows_out.append(
                        {
                            "target_normalized_symbol": sym,
                            "target_mean": float(np.mean(finite)),
                            "target_mad": raw_mad(finite),
                            "target_missing_fraction": float(1 - finite.size / tx.size),
                        }
                    )
            global_row += 1

    for sym, vecs in dup_vectors.items():
        mat = np.vstack(vecs)
        avg = np.nanmean(mat, axis=0)
        finite = avg[np.isfinite(avg)]
        rows_out.append(
            {
                "target_normalized_symbol": sym,
                "target_mean": float(np.mean(finite)),
                "target_mad": raw_mad(finite),
                "target_missing_fraction": float(1 - finite.size / avg.size),
            }
        )

    out = pd.DataFrame(rows_out)
    if out["target_normalized_symbol"].duplicated().any():
        raise RuntimeError("METABRIC marginal table still has duplicate symbols.")
    return out


def attach_matching_features(
    source: pd.DataFrame,
    target: pd.DataFrame,
    target_name: str,
) -> pd.DataFrame:
    merged = source.merge(
        target,
        left_on="Hugo_Symbol",
        right_on="target_normalized_symbol",
        how="inner",
        validate="one_to_one",
    ).copy()

    merged["source_mean_pct"] = pct_rank_average(merged["mean_log2_rsem1"])
    merged["source_mad_pct"] = pct_rank_average(merged["mad_log2_rsem1"])
    merged["target_mean_pct"] = pct_rank_average(merged["target_mean"])
    merged["target_mad_pct"] = pct_rank_average(merged["target_mad"])

    merged["source_composite"] = (
        merged["source_mean_pct"] + merged["source_mad_pct"]
    ) / 2.0
    merged["target_composite"] = (
        merged["target_mean_pct"] + merged["target_mad_pct"]
    ) / 2.0

    merged["source_composite_bin"] = quantile_bin(merged["source_composite"], q=5)
    merged["target_composite_bin"] = quantile_bin(merged["target_composite"], q=5)
    merged["matching_stratum"] = (
        merged["source_composite_bin"].astype(str)
        + "_"
        + merged["target_composite_bin"].astype(str)
    )
    merged["target"] = target_name
    return merged


def main() -> None:
    print("=" * 124)
    print("Paper 4 / TCBB - build target marginal measurability metrics for the pre-frozen mapping null")
    print("=" * 124)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Scientific guard:")
    print("  Target outcomes loaded:                              NO")
    print("  Target clinical characteristics loaded:              NO")
    print("  Target expression values loaded:                     YES")
    print("  Target gene-gene correlations calculated:            NO")
    print("  Target PCA/loadings calculated:                       NO")
    print("  Target preservation statistics calculated:           NO")
    print("  Operation: per-gene marginal mean/MAD only")
    print("=" * 124)

    for p in [
        SOURCE_METRICS,
        SOURCE_UNIVERSE,
        SCANB_EXPR,
        SCANB_PRIMARY,
        SCANB_MAP,
        METABRIC_EXPR,
        METABRIC_MAP,
        NULL_CONTRACT,
    ]:
        require(p)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    contract = json.loads(NULL_CONTRACT.read_text(encoding="utf-8"))
    if contract["status"] != "FROZEN_BEFORE_TARGET_MARGINAL_METRICS_AND_PRESERVATION":
        raise RuntimeError("Mapping-null contract is not in expected frozen state.")

    print("\n[1/4] Loading frozen source marginal metrics ...")
    source = prepare_source_features()
    print(f"  frozen source universe: {len(source):,} genes")

    print("\n[2/4] SCAN-B per-gene marginal metrics ...")
    scanb = scanb_metrics()
    print(f"  mapped SCAN-B symbols with metrics: {len(scanb):,}")

    print("\n[3/4] METABRIC per-gene marginal metrics ...")
    metabric = metabric_metrics()
    print(f"  mapped METABRIC symbols with metrics: {len(metabric):,}")

    print("\n[4/4] Building frozen 5x5 matching strata ...")
    scanb_match = attach_matching_features(source, scanb, "SCANB_GSE96058")
    metabric_match = attach_matching_features(source, metabric, "METABRIC")

    print(f"  SCAN-B frozen-10k overlap with marginal metrics:   {len(scanb_match):,}")
    print(f"  METABRIC frozen-10k overlap with marginal metrics: {len(metabric_match):,}")

    for name, df in [("SCAN-B", scanb_match), ("METABRIC", metabric_match)]:
        counts = df["matching_stratum"].value_counts()
        print(
            f"  {name} strata: n={counts.size}, "
            f"min={int(counts.min())}, median={float(counts.median()):.1f}, "
            f"max={int(counts.max())}"
        )

    scanb_out = OUT_DIR / "scanb_target_marginal_matching_features_v1.tsv"
    metabric_out = OUT_DIR / "metabric_target_marginal_matching_features_v1.tsv"
    scanb_match.to_csv(scanb_out, sep="\t", index=False)
    metabric_match.to_csv(metabric_out, sep="\t", index=False)

    result = {
        "result_id": "paper4-tcbb-target-marginal-metrics-v1",
        "script_version": SCRIPT_VERSION,
        "status": "MARGINAL_METRICS_FROZEN_NO_PRESERVATION",
        "scientific_guard": {
            "target_outcomes_loaded": False,
            "target_clinical_characteristics_loaded": False,
            "target_gene_gene_correlations_calculated": False,
            "target_pca_calculated": False,
            "target_preservation_statistics_calculated": False,
        },
        "SCANB_GSE96058": {
            "mapped_symbols_with_metrics": int(len(scanb)),
            "frozen_source_universe_overlap": int(len(scanb_match)),
        },
        "METABRIC": {
            "mapped_symbols_with_metrics": int(len(metabric)),
            "frozen_source_universe_overlap": int(len(metabric_match)),
        },
        "matching_features": contract["marginal_matching"],
    }

    json_out = OUT_DIR / "target_marginal_metrics_v1.json"
    json_out.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print("\n" + "=" * 124)
    print("03e TARGET MARGINAL METRICS: PASS")
    print("=" * 124)
    print(f"SCAN-B frozen source-universe overlap:   {len(scanb_match):,}")
    print(f"METABRIC frozen source-universe overlap: {len(metabric_match):,}")
    print()
    print("No target correlation matrix, PCA, loading, or preservation statistic was calculated.")
    print("Outputs:")
    print(f"  {scanb_out}")
    print(f"  {metabric_out}")
    print(f"  {json_out}")
    print("=" * 124)


if __name__ == "__main__":
    main()
