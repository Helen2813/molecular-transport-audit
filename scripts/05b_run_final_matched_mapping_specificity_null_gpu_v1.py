from __future__ import annotations

import gc
import hashlib
import json
import math
import re
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import torch
except Exception as exc:
    raise RuntimeError("PyTorch is required for 05b. Run the CUDA audit first.") from exc

try:
    from scipy.stats import rankdata
except Exception as exc:
    raise RuntimeError("scipy is required for exact Spearman statistics.") from exc


SCRIPT_VERSION = "05b-run-final-matched-mapping-specificity-null-gpu-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

NULL_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_mapping_null_contract_v1"
    / "structure_preserving_mapping_null_contract_v1.json"
)
PILOT_SUMMARY = (
    DATA_ROOT
    / "paper4_tcbb_mapping_null_pilot_audit_v1"
    / "mapping_null_pilot_summary_v1.tsv"
)
FINAL_EXEC_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_null_execution_v1"
    / "final_mapping_null_execution_contract_v1.json"
)
EVALUABILITY_AMENDMENT = (
    DATA_ROOT
    / "paper4_tcbb_target_statistical_evaluability_amendment_v2"
    / "target_statistical_evaluability_amendment_v2.json"
)
DIRECT_SUMMARY = (
    DATA_ROOT
    / "paper4_tcbb_primary_pooled_direct_preservation_v1"
    / "primary_pooled_direct_preservation_summary_v1.tsv"
)

SCANB_FEATURES = (
    DATA_ROOT
    / "paper4_tcbb_target_marginal_metrics_v2"
    / "scanb_target_marginal_matching_features_v2.tsv"
)
METABRIC_FEATURES = (
    DATA_ROOT
    / "paper4_tcbb_target_marginal_metrics_v2"
    / "metabric_target_marginal_matching_features_v2.tsv"
)

SCANB_EVAL = (
    DATA_ROOT
    / "paper4_tcbb_target_statistical_evaluability_amendment_v2"
    / "scanb_target_statistical_evaluability_v2.tsv"
)
METABRIC_EVAL = (
    DATA_ROOT
    / "paper4_tcbb_target_statistical_evaluability_amendment_v2"
    / "metabric_target_statistical_evaluability_v2.tsv"
)

UNIVERSE = (
    DATA_ROOT
    / "paper4_tcbb_tcga_source_universe_v1"
    / "tcga_source_gene_universe_frozen_v1.tsv"
)
WEIGHTS = (
    DATA_ROOT
    / "paper4_tcbb_frozen_source_programs_v1"
    / "tcga_frozen_source_program_weights_v1.tsv"
)

TCGA_EXPR = (
    DATA_ROOT
    / "paper4_tcbb_input_audit_v1"
    / "staged_continuous_inputs"
    / "TCGA_BRCA_PanCanAtlas2018"
    / "data_mrna_seq_v2_rsem.txt"
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
METABRIC_EXPR = (
    DATA_ROOT
    / "paper4_tcbb_input_audit_v1"
    / "staged_continuous_inputs"
    / "METABRIC"
    / "data_mrna_illumina_microarray.txt"
)

DIRECT_ROOT = DATA_ROOT / "paper4_tcbb_primary_pooled_direct_preservation_v1"
OUT_DIR = DATA_ROOT / "paper4_tcbb_final_mapping_specificity_null_v1"

VALID_PANELS_REQUIRED = 1000
MAX_PANEL_ATTEMPTS = 10_000
MAX_SPECIFICITY_EDGES = 20_000

PANEL_VALID_FRAC_LE1 = 0.90
PANEL_VALID_MAX_DISTANCE = 2

DIRECT_ALPHA = 0.05
SPECIFICITY_ALPHA = 0.05

# These bases were frozen in 03h before target preservation.
EDGE_SUBSET_BASE_SEED = 20260914
NULL_MAPPING_BASE_SEED = 20260915

# 03g pilot used the same target/program seed schedule with base 20260913.
EXPECTED_PILOT_BASE_SEED = 20260913

TARGET_INDEX = {
    "SCANB_GSE96058": 1,
    "METABRIC": 2,
}

GPU_CPU_CORR_TOL = 3e-10


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def canon_symbol(x: object) -> str:
    if x is None:
        return ""
    return " ".join(str(x).strip().split()).upper()


def module_number(program_id: str) -> int:
    m = re.search(r"M(\d+)$", str(program_id))
    if not m:
        raise RuntimeError(f"Cannot parse module number from {program_id}")
    return int(m.group(1))


def scheduled_seed(base_seed: int, target: str, program_id: str) -> int:
    """
    Mirrors the pre-result 03g pilot schedule:
      base + target_index*100000 + module_number*1000
    """
    return (
        int(base_seed)
        + TARGET_INDEX[target] * 100_000
        + module_number(program_id) * 1_000
    )


def boolish(x: object) -> bool:
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in {"1", "true", "yes", "y"}


def bh_adjust(pvalues: list[float]) -> list[float]:
    p = np.asarray(pvalues, dtype=np.float64)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    qrank = ranked * n / np.arange(1, n + 1, dtype=np.float64)
    qrank = np.minimum.accumulate(qrank[::-1])[::-1]
    qrank = np.minimum(qrank, 1.0)
    q = np.empty_like(qrank)
    q[order] = qrank
    return q.tolist()


def standardize_samples_x_genes_np(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    if not np.isfinite(x).all():
        raise RuntimeError("Non-finite values reached standardization.")
    mu = x.mean(axis=0, keepdims=True)
    sd = x.std(axis=0, ddof=1, keepdims=True)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise RuntimeError("Zero/non-finite gene SD reached standardization.")
    return (x - mu) / sd


def load_tcga_source(universe: pd.DataFrame) -> dict:
    print("\n[1/5] Replaying frozen TCGA source matrix ...")
    t0 = time.perf_counter()

    raw = pd.read_csv(TCGA_EXPR, sep="\t", low_memory=False)
    row_idx = universe["source_row_index_0based"].astype(int).to_numpy()
    selected = raw.iloc[row_idx, :]

    expected = [canon_symbol(x) for x in universe["Hugo_Symbol"]]
    observed = [canon_symbol(x) for x in selected["Hugo_Symbol"]]
    if expected != observed:
        raise RuntimeError("TCGA frozen source-universe replay mismatch.")

    values = selected.iloc[:, 2:].to_numpy(dtype=np.float64)
    if np.any(values < 0) or not np.isfinite(values).all():
        raise RuntimeError("TCGA source RSEM contains invalid values.")

    z = standardize_samples_x_genes_np(np.log2(values + 1.0).T)

    out = {
        "genes": expected,
        "gene_index": {g: i for i, g in enumerate(expected)},
        "z": z,
    }

    del raw, selected, values
    gc.collect()

    print(f"  {z.shape[0]:,} samples x {z.shape[1]:,} genes; {time.perf_counter()-t0:.1f}s")
    return out


def load_scanb_target(
    universe_set: set[str],
    primary_titles: list[str],
) -> dict:
    print("\n[2/5] Loading frozen SCAN-B primary target matrix ...")
    t0 = time.perf_counter()
    parts = []

    for chunk in pd.read_csv(
        SCANB_EXPR,
        compression="gzip",
        chunksize=1500,
        low_memory=False,
    ):
        symbol_col = chunk.columns[0]
        symbols = chunk[symbol_col].map(canon_symbol)
        mask = symbols.isin(universe_set)
        if not mask.any():
            continue

        sub = chunk.loc[mask, [symbol_col] + primary_titles]
        arr = sub[primary_titles].to_numpy(dtype=np.float64)

        fpkm = np.maximum(np.exp2(arr) - 0.1, 0.0)
        transformed = np.log2(fpkm + 1.0)

        tmp = pd.DataFrame(transformed, columns=primary_titles)
        tmp.insert(0, "Hugo_Symbol", symbols.loc[mask].to_numpy())
        parts.append(tmp)

    if not parts:
        raise RuntimeError("No SCAN-B frozen-universe genes found.")

    df = pd.concat(parts, ignore_index=True)
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    x = df.to_numpy(dtype=np.float64).T
    finite = np.isfinite(x).all(axis=0)
    sd = np.full(x.shape[1], np.nan)
    if finite.any():
        sd[finite] = x[:, finite].std(axis=0, ddof=1)
    eligible = finite & np.isfinite(sd) & (sd > 0)

    genes_all = list(df.index)
    genes = [g for g, ok in zip(genes_all, eligible) if ok]
    z = standardize_samples_x_genes_np(x[:, eligible])

    if len(genes) != 9225:
        raise RuntimeError(f"Expected 9,225 SCAN-B evaluable genes, got {len(genes)}.")

    print(f"  {z.shape[0]:,} samples x {z.shape[1]:,} evaluable genes; {time.perf_counter()-t0:.1f}s")
    return {
        "name": "SCANB_GSE96058",
        "genes": genes,
        "gene_index": {g: i for i, g in enumerate(genes)},
        "z": z,
    }


def load_metabric_target(universe_set: set[str]) -> dict:
    print("\n[2/5] Loading frozen METABRIC target matrix ...")
    t0 = time.perf_counter()
    parts = []
    sample_cols = None

    for chunk in pd.read_csv(
        METABRIC_EXPR,
        sep="\t",
        chunksize=1500,
        low_memory=False,
    ):
        if sample_cols is None:
            sample_cols = list(chunk.columns[2:])

        symbols = chunk["Hugo_Symbol"].map(canon_symbol)
        mask = symbols.isin(universe_set)
        if not mask.any():
            continue

        arr = chunk.loc[mask, sample_cols].to_numpy(dtype=np.float64)
        tmp = pd.DataFrame(arr, columns=sample_cols)
        tmp.insert(0, "Hugo_Symbol", symbols.loc[mask].to_numpy())
        parts.append(tmp)

    if not parts:
        raise RuntimeError("No METABRIC frozen-universe genes found.")

    df = pd.concat(parts, ignore_index=True)
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    x = df.to_numpy(dtype=np.float64).T
    finite = np.isfinite(x).all(axis=0)
    sd = np.full(x.shape[1], np.nan)
    if finite.any():
        sd[finite] = x[:, finite].std(axis=0, ddof=1)
    eligible = finite & np.isfinite(sd) & (sd > 0)

    genes_all = list(df.index)
    genes = [g for g, ok in zip(genes_all, eligible) if ok]
    z = standardize_samples_x_genes_np(x[:, eligible])

    if len(genes) != 8485:
        raise RuntimeError(f"Expected 8,485 METABRIC evaluable genes, got {len(genes)}.")

    print(f"  {z.shape[0]:,} samples x {z.shape[1]:,} evaluable genes; {time.perf_counter()-t0:.1f}s")
    return {
        "name": "METABRIC",
        "genes": genes,
        "gene_index": {g: i for i, g in enumerate(genes)},
        "z": z,
    }


def full_corr_gpu(z: np.ndarray, device: torch.device) -> np.ndarray:
    n = z.shape[0]
    zt = torch.as_tensor(z, dtype=torch.float64, device=device)
    ct = (zt.T @ zt) / (n - 1)
    corr = ct.detach().cpu().numpy()

    k = min(50, z.shape[1])
    cpu = (z[:, :k].T @ z[:, :k]) / (n - 1)
    diff = float(np.max(np.abs(cpu - corr[:k, :k])))
    if diff > GPU_CPU_CORR_TOL:
        raise RuntimeError(
            f"Target full-correlation CPU/GPU mismatch {diff:.3e} > {GPU_CPU_CORR_TOL:.3e}"
        )

    del zt, ct
    torch.cuda.empty_cache()
    np.fill_diagonal(corr, 1.0)

    print(f"  target full-correlation CPU↔GPU max|Δ|={diff:.3e}")
    return corr


def selected_source_correlations(
    source_z: np.ndarray,
    ii: np.ndarray,
    jj: np.ndarray,
    chunk: int = 5000,
) -> np.ndarray:
    n = source_z.shape[0]
    out = np.empty(len(ii), dtype=np.float64)

    for start in range(0, len(ii), chunk):
        end = min(start + chunk, len(ii))
        a = source_z[:, ii[start:end]]
        b = source_z[:, jj[start:end]]
        out[start:end] = np.einsum("ij,ij->j", a, b) / (n - 1)

    return out


def centered_rank(x: np.ndarray) -> tuple[np.ndarray, float]:
    r = rankdata(np.asarray(x, dtype=np.float64), method="average")
    r -= r.mean()
    ss = float(np.dot(r, r))
    if ss <= 0:
        raise RuntimeError("Degenerate rank vector.")
    return r, math.sqrt(ss)


def spearman_against_fixed_rank(
    fixed_centered_rank: np.ndarray,
    fixed_rank_norm: float,
    x: np.ndarray,
) -> float:
    r = rankdata(np.asarray(x, dtype=np.float64), method="average")
    r -= r.mean()
    norm = math.sqrt(float(np.dot(r, r)))
    if norm <= 0:
        return float("nan")
    return float(np.dot(fixed_centered_rank, r) / (fixed_rank_norm * norm))


def load_matching_features(
    target_name: str,
    eval_symbols: set[str],
) -> pd.DataFrame:
    path = SCANB_FEATURES if target_name == "SCANB_GSE96058" else METABRIC_FEATURES
    df = pd.read_csv(path, sep="\t", low_memory=False)

    required = {
        "Hugo_Symbol",
        "target_normalized_symbol",
        "source_mean_pct",
        "source_mad_pct",
        "target_mean_pct",
        "target_mad_pct",
        "source_composite_bin",
        "target_composite_bin",
    }
    if not required.issubset(df.columns):
        raise RuntimeError(
            f"{target_name}: marginal feature table missing {required - set(df.columns)}"
        )

    df = df.copy()
    df["Hugo_Symbol"] = df["Hugo_Symbol"].map(canon_symbol)
    df["target_normalized_symbol"] = df["target_normalized_symbol"].map(canon_symbol)

    # Frozen exact-symbol mapping should agree between source and target symbol fields.
    mismatch = df.loc[
        df["Hugo_Symbol"] != df["target_normalized_symbol"],
        ["Hugo_Symbol", "target_normalized_symbol"],
    ]
    if len(mismatch):
        raise RuntimeError(
            f"{target_name}: exact-symbol mapping table contains symbol mismatch; "
            f"examples={mismatch.head(10).to_dict(orient='records')}"
        )

    df = df.loc[df["Hugo_Symbol"].isin(eval_symbols)].copy()

    for c in [
        "source_mean_pct",
        "source_mad_pct",
        "target_mean_pct",
        "target_mad_pct",
    ]:
        df[c] = pd.to_numeric(df[c], errors="raise").astype(float)

    df["source_composite_bin"] = pd.to_numeric(
        df["source_composite_bin"], errors="raise"
    ).astype(int)
    df["target_composite_bin"] = pd.to_numeric(
        df["target_composite_bin"], errors="raise"
    ).astype(int)

    if df["Hugo_Symbol"].duplicated().any():
        dup = df.loc[df["Hugo_Symbol"].duplicated(False), "Hugo_Symbol"].tolist()[:20]
        raise RuntimeError(f"{target_name}: duplicate feature-table symbols: {dup}")

    if set(df["Hugo_Symbol"]) != eval_symbols:
        missing = sorted(eval_symbols - set(df["Hugo_Symbol"]))
        extra = sorted(set(df["Hugo_Symbol"]) - eval_symbols)
        raise RuntimeError(
            f"{target_name}: feature/evaluability mismatch; "
            f"missing={missing[:20]}, extra={extra[:20]}"
        )

    return df.reset_index(drop=True)


def prepare_mapping_problem(
    slot_genes: list[str],
    full_program_genes: set[str],
    features: pd.DataFrame,
    target_gene_index: dict[str, int],
) -> dict:
    feature_by_gene = features.set_index("Hugo_Symbol", drop=False)

    missing_slots = [g for g in slot_genes if g not in feature_by_gene.index]
    if missing_slots:
        raise RuntimeError(f"Missing marginal features for slot genes: {missing_slots[:20]}")

    slot = feature_by_gene.loc[slot_genes].copy().reset_index(drop=True)

    candidate = features.loc[
        ~features["Hugo_Symbol"].isin(full_program_genes)
    ].copy().reset_index(drop=True)

    missing_target_idx = [
        g for g in candidate["target_normalized_symbol"]
        if g not in target_gene_index
    ]
    if missing_target_idx:
        raise RuntimeError(
            f"Candidate genes missing from target correlation matrix: {missing_target_idx[:20]}"
        )

    candidate["target_corr_index"] = [
        target_gene_index[g] for g in candidate["target_normalized_symbol"]
    ]

    slot_bins = list(
        zip(
            slot["source_composite_bin"].astype(int),
            slot["target_composite_bin"].astype(int),
        )
    )
    cand_bins = list(
        zip(
            candidate["source_composite_bin"].astype(int),
            candidate["target_composite_bin"].astype(int),
        )
    )

    slots_by_stratum: dict[tuple[int, int], np.ndarray] = {}
    candidates_by_stratum: dict[tuple[int, int], np.ndarray] = {}

    for s in [(i, j) for i in range(1, 6) for j in range(1, 6)]:
        slots_by_stratum[s] = np.array(
            [i for i, b in enumerate(slot_bins) if b == s],
            dtype=np.int64,
        )
        candidates_by_stratum[s] = np.array(
            [i for i, b in enumerate(cand_bins) if b == s],
            dtype=np.int64,
        )

    exact_capacity_sufficient = all(
        len(candidates_by_stratum[s]) >= len(slots_by_stratum[s])
        for s in slots_by_stratum
    )

    return {
        "slot": slot,
        "candidate": candidate,
        "slot_bins": slot_bins,
        "cand_bins": cand_bins,
        "slots_by_stratum": slots_by_stratum,
        "candidates_by_stratum": candidates_by_stratum,
        "exact_capacity_sufficient": exact_capacity_sufficient,
    }


def generate_mapping(
    problem: dict,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, bool, str]:
    """
    Frozen 03d logic:
      - one-to-one unrelated target assignment;
      - same 5x5 stratum first;
      - if unavailable, expand by minimum Manhattan distance;
      - random ties;
      - no replacement.

    When every stratum has sufficient exact capacity, this fast path is
    distribution-equivalent to sequential random slot order/choice because
    strata are disjoint and no fallback can ever be triggered.
    """
    slot = problem["slot"]
    candidate = problem["candidate"]
    slots_by = problem["slots_by_stratum"]
    cand_by = problem["candidates_by_stratum"]
    m = len(slot)

    mapping_candidate_row = np.full(m, -1, dtype=np.int64)
    distances = np.full(m, -1, dtype=np.int64)

    if problem["exact_capacity_sufficient"]:
        for stratum, slot_idx in slots_by.items():
            k = len(slot_idx)
            if k == 0:
                continue

            cand_idx = cand_by[stratum]
            chosen = rng.choice(cand_idx, size=k, replace=False)
            slot_order = slot_idx[rng.permutation(k)]

            mapping_candidate_row[slot_order] = chosen
            distances[slot_order] = 0

        if np.any(mapping_candidate_row < 0):
            return mapping_candidate_row, distances, False, "exact_fast_path_incomplete"

        return mapping_candidate_row, distances, True, ""

    # General frozen fallback path for completeness.
    shuffled: dict[tuple[int, int], np.ndarray] = {}
    ptr: dict[tuple[int, int], int] = {}

    for stratum, arr in cand_by.items():
        shuffled[stratum] = arr[rng.permutation(len(arr))] if len(arr) else arr
        ptr[stratum] = 0

    slot_order = rng.permutation(m)

    for slot_i in slot_order:
        desired = problem["slot_bins"][slot_i]
        selected_stratum = None
        selected_distance = None

        for d in range(0, 9):
            tied = []
            for stratum in shuffled:
                if abs(stratum[0] - desired[0]) + abs(stratum[1] - desired[1]) != d:
                    continue
                if ptr[stratum] < len(shuffled[stratum]):
                    tied.append(stratum)

            if tied:
                selected_stratum = tied[int(rng.integers(0, len(tied)))]
                selected_distance = d
                break

        if selected_stratum is None:
            return mapping_candidate_row, distances, False, "candidate_pool_exhausted"

        p = ptr[selected_stratum]
        mapping_candidate_row[slot_i] = shuffled[selected_stratum][p]
        distances[slot_i] = int(selected_distance)
        ptr[selected_stratum] = p + 1

    return mapping_candidate_row, distances, True, ""


def mapping_quality(
    problem: dict,
    mapping_candidate_rows: np.ndarray,
    distances: np.ndarray,
) -> dict:
    slot = problem["slot"]
    candidate = problem["candidate"].iloc[mapping_candidate_rows].reset_index(drop=True)

    frac0 = float(np.mean(distances == 0))
    frac_le1 = float(np.mean(distances <= 1))
    max_dist = int(np.max(distances))
    mean_dist = float(np.mean(distances))

    return {
        "fraction_distance_0": frac0,
        "fraction_distance_le1": frac_le1,
        "maximum_distance": max_dist,
        "mean_distance": mean_dist,
        "mean_abs_source_mean_pct_diff": float(
            np.mean(
                np.abs(
                    slot["source_mean_pct"].to_numpy()
                    - candidate["source_mean_pct"].to_numpy()
                )
            )
        ),
        "mean_abs_source_mad_pct_diff": float(
            np.mean(
                np.abs(
                    slot["source_mad_pct"].to_numpy()
                    - candidate["source_mad_pct"].to_numpy()
                )
            )
        ),
        "mean_abs_target_mean_pct_diff": float(
            np.mean(
                np.abs(
                    slot["target_mean_pct"].to_numpy()
                    - candidate["target_mean_pct"].to_numpy()
                )
            )
        ),
        "mean_abs_target_mad_pct_diff": float(
            np.mean(
                np.abs(
                    slot["target_mad_pct"].to_numpy()
                    - candidate["target_mad_pct"].to_numpy()
                )
            )
        ),
        "panel_valid": bool(
            frac_le1 >= PANEL_VALID_FRAC_LE1
            and max_dist <= PANEL_VALID_MAX_DISTANCE
        ),
    }


def validate_pilot_seed_schedule(pilot: pd.DataFrame) -> None:
    bases = []

    for row in pilot.itertuples(index=False):
        target = str(row.target)
        program_id = str(row.program_id)
        seed = int(row.seed)
        base = (
            seed
            - TARGET_INDEX[target] * 100_000
            - module_number(program_id) * 1_000
        )
        bases.append(base)

    unique = sorted(set(bases))
    if unique != [EXPECTED_PILOT_BASE_SEED]:
        raise RuntimeError(
            f"03g pilot seed schedule is not the expected frozen pattern; inferred bases={unique}"
        )

    print(
        f"  03g pre-result seed schedule replay: PASS "
        f"(base={EXPECTED_PILOT_BASE_SEED})"
    )


def program_paths(target: str, program_id: str) -> dict[str, Path]:
    d = OUT_DIR / "per_program" / target
    d.mkdir(parents=True, exist_ok=True)

    return {
        "result": d / f"{program_id}_mapping_specificity_result_v1.json",
        "null": d / f"{program_id}_mapping_specificity_null_v1.npz",
        "quality": d / f"{program_id}_mapping_panel_quality_v1.tsv",
        "example": d / f"{program_id}_mapping_example_panel1_v1.tsv",
    }


def run_specificity_for_target(
    target: dict,
    target_corr: np.ndarray,
    source: dict,
    features: pd.DataFrame,
    weights: pd.DataFrame,
    direct: pd.DataFrame,
) -> list[dict]:
    target_name = target["name"]
    print("\n" + "=" * 140)
    print(f"[3/5] FINAL MATCHED-MAPPING SPECIFICITY NULL — {target_name}")
    print("=" * 140)

    results: list[dict] = []

    for p_idx, row in enumerate(
        direct.loc[direct["target"] == target_name].itertuples(index=False),
        start=1,
    ):
        program_id = str(row.program_id)
        assessable = boolish(row.primary_assessable)

        if not assessable:
            print(
                f"\n  [{p_idx:02d}/12] {program_id}: NOT ASSESSABLE — specificity not estimated"
            )
            results.append(
                {
                    "target": target_name,
                    "program_id": program_id,
                    "primary_assessable": False,
                    "specificity_estimable": False,
                    "reason": "frozen primary assessability guard failed",
                }
            )
            continue

        paths = program_paths(target_name, program_id)

        if all(paths[k].exists() for k in ["result", "null", "quality", "example"]):
            try:
                old = json.loads(paths["result"].read_text(encoding="utf-8"))
            except Exception:
                old = {}

            if (
                old.get("script_version") == SCRIPT_VERSION
                and old.get("target") == target_name
                and old.get("program_id") == program_id
                and int(old.get("valid_panels", -1)) == VALID_PANELS_REQUIRED
            ):
                print(
                    f"\n  [{p_idx:02d}/12] {program_id}: "
                    "deterministic checkpoint found — reusing completed result"
                )
                results.append(old)
                continue

        t0 = time.perf_counter()

        gene_file = (
            DIRECT_ROOT
            / "per_program"
            / target_name
            / f"{program_id}_evaluable_genes_v1.tsv"
        )
        require(gene_file)

        eg = pd.read_csv(gene_file, sep="\t", dtype=str).fillna("")
        slot_genes = [canon_symbol(x) for x in eg["Hugo_Symbol"]]
        m = len(slot_genes)

        if m != int(float(row.n_evaluable_genes)):
            raise RuntimeError(
                f"{target_name} {program_id}: evaluable gene count mismatch "
                f"{m} vs 05a {row.n_evaluable_genes}"
            )

        full_program_genes = {
            canon_symbol(x)
            for x in weights.loc[
                weights["program_id"] == program_id,
                "Hugo_Symbol",
            ]
        }

        problem = prepare_mapping_problem(
            slot_genes=slot_genes,
            full_program_genes=full_program_genes,
            features=features,
            target_gene_index=target["gene_index"],
        )

        candidate_size = len(problem["candidate"])
        expected_candidate_size = len(target["genes"]) - m
        if candidate_size != expected_candidate_size:
            raise RuntimeError(
                f"{target_name} {program_id}: candidate pool {candidate_size} "
                f"!= evaluable universe {len(target['genes'])} - slots {m}"
            )

        print(
            f"\n  [{p_idx:02d}/12] {program_id}: slots={m:,}; "
            f"candidate_pool={candidate_size:,}; "
            f"exact_capacity={'YES' if problem['exact_capacity_sufficient'] else 'NO'}"
        )

        # 03h deterministic fixed edge subset, chosen only from slot count + frozen seed.
        tri_i, tri_j = np.triu_indices(m, k=1)
        n_edges_total = len(tri_i)

        edge_seed = scheduled_seed(
            EDGE_SUBSET_BASE_SEED,
            target_name,
            program_id,
        )
        rng_edge = np.random.default_rng(edge_seed)

        if n_edges_total <= MAX_SPECIFICITY_EDGES:
            chosen = np.arange(n_edges_total, dtype=np.int64)
        else:
            chosen = np.sort(
                rng_edge.choice(
                    n_edges_total,
                    size=MAX_SPECIFICITY_EDGES,
                    replace=False,
                )
            )

        ii = tri_i[chosen]
        jj = tri_j[chosen]

        # Source structural edges for exactly the frozen evaluable source slots.
        source_cols = np.array(
            [source["gene_index"][g] for g in slot_genes],
            dtype=np.int64,
        )
        source_z = source["z"][:, source_cols]
        source_edge_subset = selected_source_correlations(source_z, ii, jj)

        # Correct target mapping on the exact same slot pairs.
        correct_target_idx = np.array(
            [target["gene_index"][g] for g in slot_genes],
            dtype=np.int64,
        )
        observed_target_edges = target_corr[
            correct_target_idx[ii],
            correct_target_idx[jj],
        ]

        source_rank, source_rank_norm = centered_rank(source_edge_subset)
        observed_rho = spearman_against_fixed_rank(
            source_rank,
            source_rank_norm,
            observed_target_edges,
        )

        full_direct_rho = float(row.rho_edge)
        if n_edges_total <= MAX_SPECIFICITY_EDGES:
            if abs(observed_rho - full_direct_rho) > 2e-12:
                raise RuntimeError(
                    f"{target_name} {program_id}: all-edge specificity rho "
                    f"{observed_rho} != 05a rho_edge {full_direct_rho}"
                )

        print(
            f"    edge subset={len(chosen):,}/{n_edges_total:,}; "
            f"observed specificity rho={observed_rho:+.6f}; "
            f"05a full rho={full_direct_rho:+.6f}"
        )

        mapping_seed = scheduled_seed(
            NULL_MAPPING_BASE_SEED,
            target_name,
            program_id,
        )
        rng = np.random.default_rng(mapping_seed)

        null_rhos: list[float] = []
        quality_rows: list[dict] = []
        attempts = 0
        example_mapping = None

        candidate_target_indices = problem["candidate"][
            "target_corr_index"
        ].to_numpy(dtype=np.int64)
        candidate_target_symbols = problem["candidate"][
            "target_normalized_symbol"
        ].to_numpy(dtype=object)

        while (
            len(null_rhos) < VALID_PANELS_REQUIRED
            and attempts < MAX_PANEL_ATTEMPTS
        ):
            attempts += 1

            mapping_rows, distances, success, failure = generate_mapping(
                problem,
                rng,
            )

            if not success:
                quality_rows.append(
                    {
                        "attempt_id": attempts,
                        "valid_panel_id": "",
                        "assignment_success": 0,
                        "failure_reason": failure,
                        "panel_valid": 0,
                    }
                )
                continue

            q = mapping_quality(
                problem,
                mapping_rows,
                distances,
            )

            mapped_target_idx = candidate_target_indices[mapping_rows]
            mapping_hash = hashlib.sha256(
                mapped_target_idx.astype("<i8", copy=False).tobytes()
            ).hexdigest()

            qrow = {
                "attempt_id": attempts,
                "valid_panel_id": (
                    len(null_rhos) + 1 if q["panel_valid"] else ""
                ),
                "assignment_success": 1,
                "failure_reason": "",
                **q,
                "mapping_sha256": mapping_hash,
            }

            if not q["panel_valid"]:
                quality_rows.append(qrow)
                continue

            target_edges = target_corr[
                mapped_target_idx[ii],
                mapped_target_idx[jj],
            ]
            rho_null = spearman_against_fixed_rank(
                source_rank,
                source_rank_norm,
                target_edges,
            )

            null_rhos.append(rho_null)
            qrow["rho_null"] = rho_null
            quality_rows.append(qrow)

            if example_mapping is None:
                candidate_rows = problem["candidate"].iloc[
                    mapping_rows
                ].reset_index(drop=True)

                example_mapping = pd.DataFrame(
                    {
                        "source_slot_index_1based": np.arange(1, m + 1),
                        "source_slot_gene": slot_genes,
                        "mapped_target_gene": candidate_rows[
                            "target_normalized_symbol"
                        ].to_numpy(),
                        "source_composite_bin_slot": problem["slot"][
                            "source_composite_bin"
                        ].to_numpy(dtype=int),
                        "target_composite_bin_slot": problem["slot"][
                            "target_composite_bin"
                        ].to_numpy(dtype=int),
                        "source_composite_bin_mapped": candidate_rows[
                            "source_composite_bin"
                        ].to_numpy(dtype=int),
                        "target_composite_bin_mapped": candidate_rows[
                            "target_composite_bin"
                        ].to_numpy(dtype=int),
                        "manhattan_distance": distances,
                    }
                )

            if len(null_rhos) % 200 == 0:
                print(
                    f"    valid panels {len(null_rhos):,}/{VALID_PANELS_REQUIRED:,} "
                    f"(attempts={attempts:,}); "
                    f"null median-so-far={np.median(null_rhos):+.5f}"
                )

        if len(null_rhos) != VALID_PANELS_REQUIRED:
            raise RuntimeError(
                f"{target_name} {program_id}: obtained {len(null_rhos)} valid panels "
                f"after {attempts} attempts; frozen requirement is {VALID_PANELS_REQUIRED} "
                f"within {MAX_PANEL_ATTEMPTS} attempts."
            )

        null = np.asarray(null_rhos, dtype=np.float64)
        extreme = int(np.sum(np.abs(null) >= abs(observed_rho)))
        p_specificity = (1 + extreme) / (VALID_PANELS_REQUIRED + 1)

        qdf = pd.DataFrame(quality_rows)
        valid_qdf = qdf.loc[qdf["panel_valid"].astype(int) == 1].copy()

        rec = {
            "script_version": SCRIPT_VERSION,
            "target": target_name,
            "program_id": program_id,
            "primary_assessable": True,
            "specificity_estimable": True,
            "n_evaluable_genes": m,
            "candidate_pool_size": candidate_size,
            "exact_stratum_capacity_sufficient": bool(
                problem["exact_capacity_sufficient"]
            ),
            "n_total_program_edges": n_edges_total,
            "n_specificity_edges": len(chosen),
            "edge_subset_seed": edge_seed,
            "mapping_seed": mapping_seed,
            "observed_specificity_rho": observed_rho,
            "direct_full_edge_rho_05a": full_direct_rho,
            "observed_minus_full_rho": observed_rho - full_direct_rho,
            "valid_panels": len(null_rhos),
            "attempts": attempts,
            "invalid_attempts": attempts - len(null_rhos),
            "mapping_specificity_empirical_p": p_specificity,
            "extreme_null_count": extreme,
            "null_mean": float(np.mean(null)),
            "null_sd": float(np.std(null, ddof=1)),
            "null_median": float(np.median(null)),
            "null_q025": float(np.quantile(null, 0.025)),
            "null_q975": float(np.quantile(null, 0.975)),
            "panel_fraction_distance_0_median": float(
                np.median(valid_qdf["fraction_distance_0"])
            ),
            "panel_fraction_distance_le1_median": float(
                np.median(valid_qdf["fraction_distance_le1"])
            ),
            "panel_maximum_distance_worst": int(
                valid_qdf["maximum_distance"].max()
            ),
            "panel_mean_distance_median": float(
                np.median(valid_qdf["mean_distance"])
            ),
            "balance_source_mean_pct_absdiff_mean": float(
                valid_qdf["mean_abs_source_mean_pct_diff"].mean()
            ),
            "balance_source_mad_pct_absdiff_mean": float(
                valid_qdf["mean_abs_source_mad_pct_diff"].mean()
            ),
            "balance_target_mean_pct_absdiff_mean": float(
                valid_qdf["mean_abs_target_mean_pct_diff"].mean()
            ),
            "balance_target_mad_pct_absdiff_mean": float(
                valid_qdf["mean_abs_target_mad_pct_diff"].mean()
            ),
            "elapsed_seconds": time.perf_counter() - t0,
        }

        np.savez_compressed(
            paths["null"],
            null_rho=null,
            edge_subset_indices_0based=chosen,
            edge_slot_i_0based=ii,
            edge_slot_j_0based=jj,
            source_edge_subset=source_edge_subset,
            correct_target_edge_subset=observed_target_edges,
        )
        qdf.to_csv(paths["quality"], sep="\t", index=False)

        if example_mapping is None:
            raise RuntimeError("No example mapping was captured.")
        example_mapping.to_csv(paths["example"], sep="\t", index=False)

        paths["result"].write_text(
            json.dumps(rec, indent=2),
            encoding="utf-8",
        )

        results.append(rec)

        print(
            f"    specificity p={p_specificity:.6f}; "
            f"null median={rec['null_median']:+.5f} "
            f"[{rec['null_q025']:+.5f},{rec['null_q975']:+.5f}]; "
            f"extreme={extreme}/{VALID_PANELS_REQUIRED}; "
            f"elapsed={rec['elapsed_seconds']:.1f}s"
        )

    return results


def apply_final_classification(
    direct: pd.DataFrame,
    specificity_results: list[dict],
) -> pd.DataFrame:
    spec_df = pd.DataFrame(specificity_results)

    estimable = spec_df.loc[
        spec_df["specificity_estimable"].map(boolish)
    ].copy()

    q_lookup: dict[tuple[str, str], float] = {}

    for target_name in sorted(estimable["target"].unique()):
        sub = estimable.loc[estimable["target"] == target_name].copy()
        p = sub["mapping_specificity_empirical_p"].astype(float).tolist()
        q = bh_adjust(p)

        for program_id, qval in zip(sub["program_id"], q):
            q_lookup[(target_name, program_id)] = qval

    rows = []

    for row in direct.itertuples(index=False):
        target = str(row.target)
        program_id = str(row.program_id)
        assessable = boolish(row.primary_assessable)

        out = {
            "target": target,
            "program_id": program_id,
            "primary_assessable": assessable,
            "n_frozen_genes": row.n_frozen_genes,
            "n_evaluable_genes": row.n_evaluable_genes,
            "coverage": row.coverage,
            "rho_edge": row.rho_edge,
            "rho_load": row.rho_load,
            "direct_edge_bh_q": row.direct_edge_bh_q,
            "direct_loading_bh_q": row.direct_loading_bh_q,
        }

        if not assessable:
            out.update(
                {
                    "mapping_specificity_rho": np.nan,
                    "mapping_specificity_p": np.nan,
                    "mapping_specificity_bh_q": np.nan,
                    "mapping_specificity_supported": 0,
                    "primary_classification": "Not assessable",
                }
            )
            rows.append(out)
            continue

        match = estimable.loc[
            (estimable["target"] == target)
            & (estimable["program_id"] == program_id)
        ]
        if len(match) != 1:
            raise RuntimeError(
                f"Expected one specificity result for {target} {program_id}, found {len(match)}"
            )

        s = match.iloc[0]
        qspec = q_lookup[(target, program_id)]

        edge_pos = boolish(row.direct_edge_positive_significant)
        load_pos = boolish(row.direct_loading_positive_significant)
        discordant = boolish(row.direct_directionally_discordant)
        spec_supported = bool(
            float(s["observed_specificity_rho"]) > 0
            and qspec < SPECIFICITY_ALPHA
        )

        if discordant:
            cls = "Directionally discordant"
        elif spec_supported and edge_pos and load_pos:
            cls = "Strong"
        elif spec_supported and (edge_pos ^ load_pos):
            cls = "Partial"
        else:
            cls = "No clear"

        out.update(
            {
                "mapping_specificity_rho": float(s["observed_specificity_rho"]),
                "mapping_specificity_p": float(
                    s["mapping_specificity_empirical_p"]
                ),
                "mapping_specificity_bh_q": qspec,
                "mapping_specificity_supported": int(spec_supported),
                "primary_classification": cls,
            }
        )
        rows.append(out)

    return pd.DataFrame(rows)


def main() -> None:
    print("=" * 142)
    print("Paper 4 / TCBB - FINAL 1,000-PANEL MATCHED-MAPPING EDGE-SPECIFICITY NULL + PRIMARY CLASSIFICATION")
    print("=" * 142)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Frozen scientific contract:")
    print("  Target outcomes/treatment loaded:                    NO")
    print("  Source program structure:                            FIXED")
    print("  Null randomizes:                                     source->target gene mapping only")
    print("  Marginal matching:                                   frozen 5x5 source/target composite strata")
    print("  Candidate pool:                                      statistically evaluable frozen target universe minus tested program")
    print(f"  Valid null panels per assessable target/program:      {VALID_PANELS_REQUIRED:,}")
    print(f"  Maximum attempts:                                    {MAX_PANEL_ATTEMPTS:,}")
    print(f"  Specificity edge subset:                             all edges if <= {MAX_SPECIFICITY_EDGES:,}, otherwise fixed {MAX_SPECIFICITY_EDGES:,}")
    print("  Loading mapping-null:                                NOT primary (frozen 03h amendment)")
    print("  BH family:                                           mapping-null p-values WITHIN target")
    print("  Final classifier:                                    frozen 04a rules")
    print("=" * 142)

    for p in [
        NULL_CONTRACT,
        PILOT_SUMMARY,
        FINAL_EXEC_CONTRACT,
        EVALUABILITY_AMENDMENT,
        DIRECT_SUMMARY,
        SCANB_FEATURES,
        METABRIC_FEATURES,
        SCANB_EVAL,
        METABRIC_EVAL,
        UNIVERSE,
        WEIGHTS,
        TCGA_EXPR,
        SCANB_EXPR,
        SCANB_PRIMARY,
        METABRIC_EXPR,
    ]:
        require(p)

    dcon = json.loads(NULL_CONTRACT.read_text(encoding="utf-8"))
    hcon = json.loads(FINAL_EXEC_CONTRACT.read_text(encoding="utf-8"))
    eva = json.loads(EVALUABILITY_AMENDMENT.read_text(encoding="utf-8"))

    if dcon.get("status") != "FROZEN_BEFORE_TARGET_MARGINAL_METRICS_AND_PRESERVATION":
        raise RuntimeError("03d mapping-null contract has unexpected status.")
    if hcon.get("status") != "FROZEN_BEFORE_FIRST_TARGET_CORRELATION_OR_PRESERVATION":
        raise RuntimeError("03h final mapping-null execution contract has unexpected status.")
    if eva.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04i evaluability amendment has unexpected status.")

    hnull = hcon["final_edge_specificity_null"]
    if int(hnull["valid_panels_required"]) != VALID_PANELS_REQUIRED:
        raise RuntimeError("05b valid-panel count differs from frozen 03h.")
    if int(hnull["maximum_panel_attempts"]) != MAX_PANEL_ATTEMPTS:
        raise RuntimeError("05b attempt cap differs from frozen 03h.")
    if int(hnull["maximum_edges_per_target_program"]) != MAX_SPECIFICITY_EDGES:
        raise RuntimeError("05b edge subset size differs from frozen 03h.")
    if int(hnull["edge_subset_seed_base"]) != EDGE_SUBSET_BASE_SEED:
        raise RuntimeError("05b edge-subset seed base differs from frozen 03h.")
    if int(hnull["null_mapping_seed_base"]) != NULL_MAPPING_BASE_SEED:
        raise RuntimeError("05b mapping seed base differs from frozen 03h.")

    pilot = pd.read_csv(PILOT_SUMMARY, sep="\t")
    validate_pilot_seed_schedule(pilot)

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable.")
    device = torch.device("cuda:0")
    prop = torch.cuda.get_device_properties(0)
    print(
        f"CUDA device: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB; "
        f"capability={prop.major}.{prop.minor}"
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    direct = pd.read_csv(DIRECT_SUMMARY, sep="\t", low_memory=False)
    if set(direct["target"]) != {"SCANB_GSE96058", "METABRIC"}:
        raise RuntimeError(f"Unexpected 05a targets: {sorted(set(direct['target']))}")

    universe = pd.read_csv(UNIVERSE, sep="\t", dtype=str).fillna("")
    universe["source_gene_rank_by_MAD"] = pd.to_numeric(
        universe["source_gene_rank_by_MAD"], errors="raise"
    ).astype(int)
    universe = universe.sort_values("source_gene_rank_by_MAD").reset_index(drop=True)

    if len(universe) != 10_000:
        raise RuntimeError(f"Expected 10,000 source-universe genes, got {len(universe)}.")

    universe_genes = [canon_symbol(x) for x in universe["Hugo_Symbol"]]
    universe_set = set(universe_genes)
    if len(universe_set) != 10_000:
        raise RuntimeError("Source universe is not unique after canonicalization.")

    weights = pd.read_csv(WEIGHTS, sep="\t", dtype=str).fillna("")
    weights["Hugo_Symbol"] = weights["Hugo_Symbol"].map(canon_symbol)

    primary = pd.read_csv(SCANB_PRIMARY, sep="\t", dtype=str).fillna("")
    primary_titles = list(primary["primary_title"])
    if len(primary_titles) != 3273 or len(set(primary_titles)) != 3273:
        raise RuntimeError("SCAN-B frozen primary manifest is not exactly 3,273 unique profiles.")

    scanb_eval = pd.read_csv(SCANB_EVAL, sep="\t", dtype=str).fillna("")
    met_eval = pd.read_csv(METABRIC_EVAL, sep="\t", dtype=str).fillna("")

    scanb_eval_symbols = {
        canon_symbol(x)
        for x in scanb_eval.loc[
            scanb_eval["statistically_evaluable"].astype(str) == "1",
            "Hugo_Symbol",
        ]
    }
    met_eval_symbols = {
        canon_symbol(x)
        for x in met_eval.loc[
            met_eval["statistically_evaluable"].astype(str) == "1",
            "Hugo_Symbol",
        ]
    }

    if len(scanb_eval_symbols) != 9225:
        raise RuntimeError(f"Expected 9,225 SCAN-B evaluable feature genes, got {len(scanb_eval_symbols)}.")
    if len(met_eval_symbols) != 8485:
        raise RuntimeError(f"Expected 8,485 METABRIC evaluable feature genes, got {len(met_eval_symbols)}.")

    scanb_features = load_matching_features(
        "SCANB_GSE96058",
        scanb_eval_symbols,
    )
    met_features = load_matching_features(
        "METABRIC",
        met_eval_symbols,
    )

    source = load_tcga_source(universe)

    all_spec_results: list[dict] = []

    # Process one target at a time to keep the full 8-9k gene correlation matrix memory bounded.
    for target_name in ["SCANB_GSE96058", "METABRIC"]:
        target = (
            load_scanb_target(universe_set, primary_titles)
            if target_name == "SCANB_GSE96058"
            else load_metabric_target(universe_set)
        )
        features = scanb_features if target_name == "SCANB_GSE96058" else met_features

        print(
            f"\n  Computing one full {len(target['genes']):,}x{len(target['genes']):,} "
            f"target Pearson matrix on GPU for {target_name} ..."
        )
        t_corr = time.perf_counter()
        target_corr = full_corr_gpu(target["z"], device)
        print(
            f"  full target correlation matrix ready in "
            f"{time.perf_counter()-t_corr:.1f}s"
        )

        all_spec_results.extend(
            run_specificity_for_target(
                target=target,
                target_corr=target_corr,
                source=source,
                features=features,
                weights=weights,
                direct=direct,
            )
        )

        del target_corr, target
        torch.cuda.empty_cache()
        gc.collect()

    print("\n[4/5] Applying frozen within-target BH correction to mapping-specificity p-values ...")
    classification = apply_final_classification(
        direct=direct,
        specificity_results=all_spec_results,
    )

    summary_path = OUT_DIR / "primary_pooled_final_classification_v1.tsv"
    classification.to_csv(summary_path, sep="\t", index=False)

    spec_df = pd.DataFrame(all_spec_results)
    spec_path = OUT_DIR / "mapping_specificity_summary_v1.tsv"
    spec_df.to_csv(spec_path, sep="\t", index=False)

    print("\n[5/5] Final frozen primary classifications ...")
    class_counts = {}

    for target_name in ["SCANB_GSE96058", "METABRIC"]:
        sub = classification.loc[classification["target"] == target_name].copy()
        counts = sub["primary_classification"].value_counts().to_dict()
        class_counts[target_name] = {str(k): int(v) for k, v in counts.items()}

        print(f"\n{target_name}:")
        for r in sub.itertuples(index=False):
            if r.primary_classification == "Not assessable":
                print(
                    f"  {r.program_id}: NOT ASSESSABLE "
                    f"(eligible={int(float(r.n_evaluable_genes))}/{int(float(r.n_frozen_genes))})"
                )
            else:
                print(
                    f"  {r.program_id}: "
                    f"rho_edge={float(r.rho_edge):+.4f}; "
                    f"rho_load={float(r.rho_load):+.4f}; "
                    f"rho_specificity={float(r.mapping_specificity_rho):+.4f}; "
                    f"q_specificity={float(r.mapping_specificity_bh_q):.6f}; "
                    f"CLASS={r.primary_classification}"
                )

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "FINAL_PRIMARY_POOLED_CLASSIFICATION_COMPLETE",
        "valid_panels_per_assessable_combination": VALID_PANELS_REQUIRED,
        "maximum_specificity_edges": MAX_SPECIFICITY_EDGES,
        "edge_subset_seed_base": EDGE_SUBSET_BASE_SEED,
        "null_mapping_seed_base": NULL_MAPPING_BASE_SEED,
        "seed_derivation": (
            "base + target_index*100000 + module_number*1000; "
            "same schedule form as pre-result 03g pilot"
        ),
        "classification_counts": class_counts,
        "classification_file": str(summary_path),
        "specificity_summary_file": str(spec_path),
        "target_outcomes_loaded": False,
        "loading_mapping_null_primary": False,
    }

    master_path = OUT_DIR / "final_mapping_specificity_and_classification_v1.json"
    master_path.write_text(json.dumps(master, indent=2), encoding="utf-8")

    print("\n" + "=" * 142)
    print("05b FINAL MATCHED-MAPPING SPECIFICITY NULL: COMPLETE")
    print("=" * 142)
    print("The frozen primary pooled Strong/Partial/No-clear classifications are now available.")
    print("No target outcomes or treatment variables were loaded.")
    print()
    print(f"Classification TSV: {summary_path}")
    print(f"Specificity TSV:    {spec_path}")
    print(f"Master JSON:        {master_path}")
    print("=" * 142)


if __name__ == "__main__":
    main()
