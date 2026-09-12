from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_VERSION = "03e-build-target-marginal-metrics-v2-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

# Authoritative frozen source universe. This table already contains the exact
# source-row choice after duplicate-symbol resolution plus the source marginal
# metrics used for matching.
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

# New output directory: v1 failed before any target expression was read.
OUT_DIR = DATA_ROOT / "paper4_tcbb_target_marginal_metrics_v2"

CHUNK_ROWS = 128


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def normalize_symbol(x: object) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def raw_mad(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    med = np.median(x)
    return float(np.median(np.abs(x - med)))


def pct_rank_average(x: pd.Series) -> pd.Series:
    return x.rank(method="average", pct=True)


def quantile_bin(score: pd.Series, q: int = 5) -> pd.Series:
    # Rank first so binning is deterministic even with tied composite values.
    ranks = score.rank(method="first")
    return pd.qcut(ranks, q=q, labels=False, duplicates="raise").astype(int) + 1


def prepare_source_features() -> pd.DataFrame:
    """
    Use the authoritative frozen 10k source-universe manifest directly.

    v1 incorrectly attempted a one-to-one merge back to the pre-deduplication
    20,531-row audit table, which legitimately contains duplicate Hugo symbols.
    The frozen universe already records the exact retained row and its source
    mean/MAD, so no merge or re-deduplication is scientifically appropriate.
    """
    universe = pd.read_csv(
        SOURCE_UNIVERSE,
        sep="\t",
        dtype={"Hugo_Symbol": str, "Entrez_Gene_Id": str},
        low_memory=False,
    )

    required_cols = {
        "Hugo_Symbol",
        "source_row_index_0based",
        "mean_log2_rsem1",
        "mad_log2_rsem1",
    }
    missing = required_cols.difference(universe.columns)
    if missing:
        raise RuntimeError(
            f"Frozen source universe missing required columns: {sorted(missing)}"
        )

    if len(universe) != 10000:
        raise RuntimeError(
            f"Expected exactly 10,000 frozen source genes; found {len(universe)}"
        )

    out = universe[
        [
            "Hugo_Symbol",
            "source_row_index_0based",
            "mean_log2_rsem1",
            "mad_log2_rsem1",
        ]
    ].copy()

    out["Hugo_Symbol"] = out["Hugo_Symbol"].map(normalize_symbol)

    if out["Hugo_Symbol"].eq("").any():
        raise RuntimeError("Blank symbol in frozen source universe.")
    if out["Hugo_Symbol"].duplicated().any():
        dup = out.loc[
            out["Hugo_Symbol"].duplicated(keep=False), "Hugo_Symbol"
        ].tolist()
        raise RuntimeError(
            f"Frozen source universe is not unique after normalization: {dup[:20]}"
        )

    for col in ["mean_log2_rsem1", "mad_log2_rsem1"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
        if out[col].isna().any():
            raise RuntimeError(f"Non-numeric/missing frozen source metric: {col}")

    return out


def parse_mapping_rows(mapping_path: Path) -> tuple[dict[int, str], set[str]]:
    mapping = pd.read_csv(mapping_path, sep="\t", dtype=str)
    row_to_symbol: dict[int, str] = {}
    duplicate_symbols: set[str] = set()

    for _, r in mapping.iterrows():
        sym = normalize_symbol(r["target_normalized_symbol"])
        idxs = [int(x) for x in str(r["target_row_indices_0based"]).split(",")]
        for idx in idxs:
            if idx in row_to_symbol:
                raise RuntimeError(f"Target row index mapped twice: {idx}")
            row_to_symbol[idx] = sym
        if len(idxs) > 1:
            duplicate_symbols.add(sym)

    return row_to_symbol, duplicate_symbols


def scanb_metrics() -> pd.DataFrame:
    primary = pd.read_csv(SCANB_PRIMARY, sep="\t", dtype=str)
    primary_cols = primary["primary_title"].tolist()

    if len(primary_cols) != 3273:
        raise RuntimeError(
            f"Expected 3,273 frozen primary SCAN-B profiles; found {len(primary_cols)}"
        )
    if len(primary_cols) != len(set(primary_cols)):
        raise RuntimeError("Duplicate frozen SCAN-B primary profile title.")

    row_to_symbol, dup_symbols = parse_mapping_rows(SCANB_MAP)

    header = pd.read_csv(SCANB_EXPR, compression="gzip", nrows=0)
    all_cols = list(header.columns)
    id_col = all_cols[0]

    missing = [c for c in primary_cols if c not in all_cols]
    if missing:
        raise RuntimeError(f"SCAN-B primary columns missing: {missing[:10]}")

    rows_out: list[dict] = []
    dup_vectors: dict[str, list[np.ndarray]] = defaultdict(list)
    global_row = 0

    usecols = [id_col] + primary_cols
    for chunk_i, chunk in enumerate(
        pd.read_csv(
            SCANB_EXPR,
            compression="gzip",
            usecols=usecols,
            chunksize=CHUNK_ROWS,
            low_memory=False,
        ),
        start=1,
    ):
        values = chunk[primary_cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)

        for local_i in range(len(chunk)):
            if global_row in row_to_symbol:
                sym = row_to_symbol[global_row]
                y = values[local_i, :]

                # Frozen SCAN-B primary representation:
                # published y = log2(FPKM + 0.1)
                # -> FPKM -> log2(FPKM + 1)
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
                            "target_missing_fraction": float(
                                1 - finite.size / tx.size
                            ),
                        }
                    )
            global_row += 1

        if chunk_i % 50 == 0:
            print(f"    SCAN-B processed {global_row:,} gene rows ...")

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
    row_to_symbol, dup_symbols = parse_mapping_rows(METABRIC_MAP)

    header = pd.read_csv(METABRIC_EXPR, sep="\t", nrows=0)
    all_cols = list(header.columns)
    sample_cols = all_cols[2:]

    if len(sample_cols) != 1980:
        raise RuntimeError(
            f"Expected 1,980 METABRIC profiles; found {len(sample_cols)}"
        )

    rows_out: list[dict] = []
    dup_vectors: dict[str, list[np.ndarray]] = defaultdict(list)
    global_row = 0

    for chunk_i, chunk in enumerate(
        pd.read_csv(
            METABRIC_EXPR,
            sep="\t",
            chunksize=CHUNK_ROWS,
            low_memory=False,
        ),
        start=1,
    ):
        values = chunk[sample_cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)

        for local_i in range(len(chunk)):
            if global_row in row_to_symbol:
                sym = row_to_symbol[global_row]
                tx = values[local_i, :]

                # Frozen METABRIC primary representation is as provided:
                # continuous log2 Illumina HT-12 v3 intensity.
                if sym in dup_symbols:
                    dup_vectors[sym].append(tx)
                else:
                    finite = tx[np.isfinite(tx)]
                    rows_out.append(
                        {
                            "target_normalized_symbol": sym,
                            "target_mean": float(np.mean(finite)),
                            "target_mad": raw_mad(finite),
                            "target_missing_fraction": float(
                                1 - finite.size / tx.size
                            ),
                        }
                    )
            global_row += 1

        if chunk_i % 50 == 0:
            print(f"    METABRIC processed {global_row:,} gene rows ...")

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

    merged["source_composite_bin"] = quantile_bin(
        merged["source_composite"], q=5
    )
    merged["target_composite_bin"] = quantile_bin(
        merged["target_composite"], q=5
    )
    merged["matching_stratum"] = (
        merged["source_composite_bin"].astype(str)
        + "_"
        + merged["target_composite_bin"].astype(str)
    )
    merged["target"] = target_name
    return merged


def main() -> None:
    print("=" * 124)
    print(
        "Paper 4 / TCBB - build target marginal measurability metrics "
        "for the pre-frozen mapping null"
    )
    print("=" * 124)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("v2 correction:")
    print(
        "  v1 failed before opening target expression because it attempted "
        "a one-to-one merge against the pre-dedup source audit table."
    )
    print(
        "  v2 uses the authoritative frozen 10,000-gene source-universe "
        "manifest directly; no scientific rule is changed."
    )
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
    if (
        contract["status"]
        != "FROZEN_BEFORE_TARGET_MARGINAL_METRICS_AND_PRESERVATION"
    ):
        raise RuntimeError(
            "Mapping-null contract is not in expected frozen state."
        )

    print("\n[1/4] Loading authoritative frozen source marginal metrics ...")
    source = prepare_source_features()
    print(f"  frozen source universe: {len(source):,} unique genes")

    print("\n[2/4] SCAN-B per-gene marginal metrics ...")
    scanb = scanb_metrics()
    print(f"  mapped SCAN-B symbols with metrics: {len(scanb):,}")

    print("\n[3/4] METABRIC per-gene marginal metrics ...")
    metabric = metabric_metrics()
    print(f"  mapped METABRIC symbols with metrics: {len(metabric):,}")

    print("\n[4/4] Building frozen 5x5 matching strata ...")
    scanb_match = attach_matching_features(
        source, scanb, "SCANB_GSE96058"
    )
    metabric_match = attach_matching_features(
        source, metabric, "METABRIC"
    )

    # These counts must replay the already-frozen identifier overlaps from 03a.
    if len(scanb_match) != 9225:
        raise RuntimeError(
            f"Expected SCAN-B frozen-universe overlap 9,225; found {len(scanb_match)}"
        )
    if len(metabric_match) != 8490:
        raise RuntimeError(
            f"Expected METABRIC frozen-universe overlap 8,490; found {len(metabric_match)}"
        )

    print(
        f"  SCAN-B frozen-10k overlap with marginal metrics:   {len(scanb_match):,}"
    )
    print(
        f"  METABRIC frozen-10k overlap with marginal metrics: "
        f"{len(metabric_match):,}"
    )

    stratum_summaries = {}
    for name, df in [
        ("SCAN-B", scanb_match),
        ("METABRIC", metabric_match),
    ]:
        counts = df["matching_stratum"].value_counts()
        stratum_summaries[name] = {
            "n_strata": int(counts.size),
            "min": int(counts.min()),
            "median": float(counts.median()),
            "max": int(counts.max()),
        }
        print(
            f"  {name} strata: n={counts.size}, "
            f"min={int(counts.min())}, "
            f"median={float(counts.median()):.1f}, "
            f"max={int(counts.max())}"
        )

    scanb_out = (
        OUT_DIR / "scanb_target_marginal_matching_features_v2.tsv"
    )
    metabric_out = (
        OUT_DIR / "metabric_target_marginal_matching_features_v2.tsv"
    )
    scanb_match.to_csv(scanb_out, sep="\t", index=False)
    metabric_match.to_csv(metabric_out, sep="\t", index=False)

    result = {
        "result_id": "paper4-tcbb-target-marginal-metrics-v2",
        "script_version": SCRIPT_VERSION,
        "supersedes_failed_script": (
            "03e-build-target-marginal-metrics-v1-no-cli"
        ),
        "v1_failure_stage": (
            "source-only pre-target merge; no target expression had been opened"
        ),
        "status": "MARGINAL_METRICS_FROZEN_NO_PRESERVATION",
        "scientific_guard": {
            "target_outcomes_loaded": False,
            "target_clinical_characteristics_loaded": False,
            "target_gene_gene_correlations_calculated": False,
            "target_pca_calculated": False,
            "target_preservation_statistics_calculated": False,
        },
        "source_feature_authority": (
            "frozen 10,000-gene source-universe manifest after source-only "
            "duplicate-symbol resolution"
        ),
        "SCANB_GSE96058": {
            "mapped_symbols_with_metrics": int(len(scanb)),
            "frozen_source_universe_overlap": int(len(scanb_match)),
        },
        "METABRIC": {
            "mapped_symbols_with_metrics": int(len(metabric)),
            "frozen_source_universe_overlap": int(len(metabric_match)),
        },
        "strata": stratum_summaries,
        "matching_features": contract["marginal_matching"],
    }

    json_out = OUT_DIR / "target_marginal_metrics_v2.json"
    json_out.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print("\n" + "=" * 124)
    print("03e v2 TARGET MARGINAL METRICS: PASS")
    print("=" * 124)
    print(
        f"SCAN-B frozen source-universe overlap:   {len(scanb_match):,}"
    )
    print(
        f"METABRIC frozen source-universe overlap: {len(metabric_match):,}"
    )
    print()
    print(
        "No target correlation matrix, PCA, loading, or preservation "
        "statistic was calculated."
    )
    print("Outputs:")
    print(f"  {scanb_out}")
    print(f"  {metabric_out}")
    print(f"  {json_out}")
    print("=" * 124)


if __name__ == "__main__":
    main()
