from __future__ import annotations

import gc
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import torch
except Exception as exc:
    raise RuntimeError("PyTorch is required for 05d.") from exc

try:
    from scipy.stats import rankdata
except Exception as exc:
    raise RuntimeError("scipy is required for Spearman concordance.") from exc


SCRIPT_VERSION = "05d-run-controlled-mapping-degradation-gpu-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

EXEC_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_controlled_degradation_execution_contract_v1"
    / "controlled_degradation_execution_contract_v1.json"
)
OPERATING_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_preservation_operating_contract_v1"
    / "preservation_operating_characteristics_contract_v1.json"
)
FINAL_CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
    / "primary_pooled_final_classification_v2.tsv"
)
SPECIFICITY_SUMMARY = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
    / "mapping_specificity_summary_v2.tsv"
)
SPECIFICITY_ROOT = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
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
    / "paper4_tcbb_exact_variation_evaluability_correction_v1"
    / "scanb_target_exact_variation_evaluability_v1.tsv"
)
METABRIC_EVAL = (
    DATA_ROOT
    / "paper4_tcbb_exact_variation_evaluability_correction_v1"
    / "metabric_target_exact_variation_evaluability_v1.tsv"
)

WEIGHTS = (
    DATA_ROOT
    / "paper4_tcbb_frozen_source_programs_v1"
    / "tcga_frozen_source_program_weights_v1.tsv"
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

DIRECT_ROOT = DATA_ROOT / "paper4_tcbb_primary_pooled_direct_preservation_v2"
OUT_DIR = DATA_ROOT / "paper4_tcbb_controlled_mapping_degradation_v1"

FRACTIONS = [0.0, 0.10, 0.25, 0.50, 0.75, 1.0]
REPLICATES = 100
MAX_ATTEMPTS_PER_FRACTION = 1000
BASE_SEED = 20260916

MIN_FRAC_LE1 = 0.90
MAX_DISTANCE = 2
TREND_GATE = -0.80
FALSE_SPECIFICITY_MAX = 0.10
ALPHA = 0.05
RHO_REPLAY_TOL = 5e-12

TARGET_INDEX = {
    "SCANB_GSE96058": 1,
    "METABRIC": 2,
}


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
        raise RuntimeError(f"Cannot parse module number from {program_id}.")
    return int(m.group(1))


def boolish(x: object) -> bool:
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in {"1", "true", "yes", "y"}


def centered_rank(x: np.ndarray) -> tuple[np.ndarray, float]:
    r = rankdata(np.asarray(x, dtype=np.float64), method="average")
    r -= r.mean()
    ss = float(np.dot(r, r))
    if ss <= 0:
        raise RuntimeError("Degenerate rank vector.")
    return r, math.sqrt(ss)


def spearman_against_fixed(
    fixed_rank: np.ndarray,
    fixed_norm: float,
    x: np.ndarray,
) -> float:
    r = rankdata(np.asarray(x, dtype=np.float64), method="average")
    r -= r.mean()
    norm = math.sqrt(float(np.dot(r, r)))
    if norm <= 0:
        return float("nan")
    return float(np.dot(fixed_rank, r) / (fixed_norm * norm))


def standardize_samples_x_genes_np(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    if not np.isfinite(x).all():
        raise RuntimeError("Non-finite values reached standardization.")
    mu = x.mean(axis=0, keepdims=True)
    sd = x.std(axis=0, ddof=1, keepdims=True)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise RuntimeError("Zero/non-finite SD reached standardization.")
    return (x - mu) / sd


def load_scanb_target(
    universe_set: set[str],
    primary_titles: list[str],
) -> dict:
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
        raise RuntimeError("No SCAN-B corrected-universe genes found.")

    df = pd.concat(parts, ignore_index=True)
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    x = df.to_numpy(dtype=np.float64).T
    finite = np.isfinite(x).all(axis=0)
    exact_min = np.full(x.shape[1], np.nan)
    exact_max = np.full(x.shape[1], np.nan)
    exact_min[finite] = np.min(x[:, finite], axis=0)
    exact_max[finite] = np.max(x[:, finite], axis=0)
    eligible = finite & (exact_max > exact_min)

    genes_all = list(df.index)
    genes = [g for g, ok in zip(genes_all, eligible) if ok]
    z = standardize_samples_x_genes_np(x[:, eligible])

    if len(genes) != 9220:
        raise RuntimeError(f"Expected 9,220 corrected SCAN-B genes, got {len(genes)}.")

    return {
        "name": "SCANB_GSE96058",
        "genes": genes,
        "gene_index": {g: i for i, g in enumerate(genes)},
        "z": z,
    }


def load_metabric_target(universe_set: set[str]) -> dict:
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
        raise RuntimeError("No METABRIC corrected-universe genes found.")

    df = pd.concat(parts, ignore_index=True)
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    x = df.to_numpy(dtype=np.float64).T
    finite = np.isfinite(x).all(axis=0)
    exact_min = np.full(x.shape[1], np.nan)
    exact_max = np.full(x.shape[1], np.nan)
    exact_min[finite] = np.min(x[:, finite], axis=0)
    exact_max[finite] = np.max(x[:, finite], axis=0)
    eligible = finite & (exact_max > exact_min)

    genes_all = list(df.index)
    genes = [g for g, ok in zip(genes_all, eligible) if ok]
    z = standardize_samples_x_genes_np(x[:, eligible])

    if len(genes) != 8485:
        raise RuntimeError(f"Expected 8,485 corrected METABRIC genes, got {len(genes)}.")

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
    if diff > 3e-10:
        raise RuntimeError(f"CPU/GPU target correlation mismatch: {diff:.3e}")

    del zt, ct
    torch.cuda.empty_cache()
    np.fill_diagonal(corr, 1.0)
    print(f"  target full-correlation CPU↔GPU max|Δ|={diff:.3e}")
    return corr


def load_matching_features(
    target_name: str,
    eval_symbols: set[str],
) -> pd.DataFrame:
    path = SCANB_FEATURES if target_name == "SCANB_GSE96058" else METABRIC_FEATURES
    df = pd.read_csv(path, sep="\t", low_memory=False)

    df["Hugo_Symbol"] = df["Hugo_Symbol"].map(canon_symbol)
    df["target_normalized_symbol"] = df["target_normalized_symbol"].map(canon_symbol)

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
        raise RuntimeError(f"{target_name}: duplicate marginal-feature symbols.")

    if set(df["Hugo_Symbol"]) != eval_symbols:
        missing = sorted(eval_symbols - set(df["Hugo_Symbol"]))
        raise RuntimeError(
            f"{target_name}: corrected evaluability/marginal-feature mismatch; "
            f"missing={missing[:20]}"
        )

    if not (df["Hugo_Symbol"] == df["target_normalized_symbol"]).all():
        raise RuntimeError(f"{target_name}: exact-symbol mismatch in marginal features.")

    return df.reset_index(drop=True)


def prepare_problem(
    slot_genes: list[str],
    corrupted_slot_idx: np.ndarray,
    full_program_genes: set[str],
    features: pd.DataFrame,
    target_gene_index: dict[str, int],
) -> dict:
    by_gene = features.set_index("Hugo_Symbol", drop=False)
    selected_genes = [slot_genes[i] for i in corrupted_slot_idx]
    slot = by_gene.loc[selected_genes].copy().reset_index(drop=True)

    candidate = features.loc[
        ~features["Hugo_Symbol"].isin(full_program_genes)
    ].copy().reset_index(drop=True)

    candidate["target_corr_index"] = [
        target_gene_index[g]
        for g in candidate["target_normalized_symbol"]
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

    strata = [(i, j) for i in range(1, 6) for j in range(1, 6)]
    slots_by = {
        s: np.array([i for i, b in enumerate(slot_bins) if b == s], dtype=np.int64)
        for s in strata
    }
    cand_by = {
        s: np.array([i for i, b in enumerate(cand_bins) if b == s], dtype=np.int64)
        for s in strata
    }
    exact_capacity = all(
        len(cand_by[s]) >= len(slots_by[s])
        for s in strata
    )

    return {
        "slot": slot,
        "candidate": candidate,
        "slot_bins": slot_bins,
        "slots_by": slots_by,
        "cand_by": cand_by,
        "exact_capacity": exact_capacity,
    }


def generate_mapping(
    problem: dict,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, bool, str]:
    m = len(problem["slot"])
    mapping_rows = np.full(m, -1, dtype=np.int64)
    distances = np.full(m, -1, dtype=np.int64)

    if problem["exact_capacity"]:
        for stratum, slot_idx in problem["slots_by"].items():
            k = len(slot_idx)
            if k == 0:
                continue
            cand_idx = problem["cand_by"][stratum]
            chosen = rng.choice(cand_idx, size=k, replace=False)
            slot_order = slot_idx[rng.permutation(k)]
            mapping_rows[slot_order] = chosen
            distances[slot_order] = 0

        if np.any(mapping_rows < 0):
            return mapping_rows, distances, False, "exact_fast_path_incomplete"
        return mapping_rows, distances, True, ""

    shuffled = {}
    ptr = {}
    for stratum, arr in problem["cand_by"].items():
        shuffled[stratum] = arr[rng.permutation(len(arr))] if len(arr) else arr
        ptr[stratum] = 0

    order = rng.permutation(m)

    for slot_i in order:
        desired = problem["slot_bins"][slot_i]
        selected = None
        selected_distance = None

        for d in range(0, 9):
            tied = []
            for stratum in shuffled:
                if abs(stratum[0] - desired[0]) + abs(stratum[1] - desired[1]) != d:
                    continue
                if ptr[stratum] < len(shuffled[stratum]):
                    tied.append(stratum)

            if tied:
                selected = tied[int(rng.integers(0, len(tied)))]
                selected_distance = d
                break

        if selected is None:
            return mapping_rows, distances, False, "candidate_pool_exhausted"

        p = ptr[selected]
        mapping_rows[slot_i] = shuffled[selected][p]
        distances[slot_i] = int(selected_distance)
        ptr[selected] = p + 1

    return mapping_rows, distances, True, ""


def quality(distances: np.ndarray) -> dict:
    frac0 = float(np.mean(distances == 0))
    frac_le1 = float(np.mean(distances <= 1))
    maxd = int(np.max(distances))
    meand = float(np.mean(distances))

    return {
        "fraction_distance_0": frac0,
        "fraction_distance_le1": frac_le1,
        "maximum_distance": maxd,
        "mean_distance": meand,
        "panel_valid": bool(
            frac_le1 >= MIN_FRAC_LE1
            and maxd <= MAX_DISTANCE
        ),
    }


def rng_for_attempt(
    target: str,
    program_id: str,
    fraction: float,
    attempt_id: int,
) -> np.random.Generator:
    fraction_code = int(round(1000.0 * fraction))
    ss = np.random.SeedSequence(
        [
            BASE_SEED,
            TARGET_INDEX[target],
            module_number(program_id),
            fraction_code,
            int(attempt_id),
        ]
    )
    return np.random.default_rng(ss)


def round_half_up(x: float) -> int:
    return int(math.floor(x + 0.5))


def run_combo(
    target: dict,
    target_corr: np.ndarray,
    features: pd.DataFrame,
    weights: pd.DataFrame,
    class_row: dict,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    target_name = str(class_row["target"])
    program_id = str(class_row["program_id"])

    gene_file = (
        DIRECT_ROOT
        / "per_program"
        / target_name
        / f"{program_id}_evaluable_genes_v2.tsv"
    )
    require(gene_file)
    eg = pd.read_csv(gene_file, sep="\t", dtype=str).fillna("")
    slot_genes = [canon_symbol(x) for x in eg["Hugo_Symbol"]]
    m = len(slot_genes)

    full_program_genes = {
        canon_symbol(x)
        for x in weights.loc[
            weights["program_id"] == program_id,
            "Hugo_Symbol",
        ]
    }

    null_path = (
        SPECIFICITY_ROOT
        / "per_program"
        / target_name
        / f"{program_id}_mapping_specificity_null_v2.npz"
    )
    require(null_path)
    nz = np.load(null_path)

    frozen_null = np.asarray(nz["null_rho"], dtype=np.float64)
    ii = np.asarray(nz["edge_slot_i_0based"], dtype=np.int64)
    jj = np.asarray(nz["edge_slot_j_0based"], dtype=np.int64)
    source_edges = np.asarray(nz["source_edge_subset"], dtype=np.float64)

    if len(ii) != len(jj) or len(ii) != len(source_edges):
        raise RuntimeError(f"{target_name} {program_id}: frozen subset length mismatch.")
    if int(max(ii.max(), jj.max())) >= m:
        raise RuntimeError(f"{target_name} {program_id}: frozen edge subset exceeds slot count.")

    source_rank, source_norm = centered_rank(source_edges)

    correct_target_idx = np.array(
        [target["gene_index"][g] for g in slot_genes],
        dtype=np.int64,
    )
    correct_edges = target_corr[
        correct_target_idx[ii],
        correct_target_idx[jj],
    ]
    rho0 = spearman_against_fixed(
        source_rank,
        source_norm,
        correct_edges,
    )

    expected_rho0 = float(class_row["mapping_specificity_rho"])
    if abs(rho0 - expected_rho0) > RHO_REPLAY_TOL:
        raise RuntimeError(
            f"{target_name} {program_id}: zero-fraction rho replay mismatch "
            f"{rho0} vs corrected 05b v2 {expected_rho0}"
        )

    null_q025 = float(np.quantile(frozen_null, 0.025))
    null_q975 = float(np.quantile(frozen_null, 0.975))

    replicate_rows = [
        {
            "target": target_name,
            "program_id": program_id,
            "fraction": 0.0,
            "replicate_id": 0,
            "attempt_id": 0,
            "n_slots": m,
            "n_corrupted_slots": 0,
            "realized_corruption_fraction": 0.0,
            "rho_edge": rho0,
            "panel_valid": 1,
            "fraction_distance_0": np.nan,
            "fraction_distance_le1": np.nan,
            "maximum_distance": np.nan,
            "mean_distance": np.nan,
            "nominal_mapping_null_p": np.nan,
        }
    ]
    fraction_summaries = [
        {
            "target": target_name,
            "program_id": program_id,
            "fraction": 0.0,
            "n_replicates": 1,
            "n_corrupted_slots": 0,
            "realized_corruption_fraction": 0.0,
            "rho_median": rho0,
            "rho_q025": rho0,
            "rho_q975": rho0,
            "attempts": 0,
            "invalid_attempts": 0,
        }
    ]

    for fraction in FRACTIONS[1:]:
        k = max(1, min(m, round_half_up(fraction * m)))
        valid_rhos = []
        valid_rows = []
        attempts = 0

        while (
            len(valid_rhos) < REPLICATES
            and attempts < MAX_ATTEMPTS_PER_FRACTION
        ):
            attempts += 1
            rng = rng_for_attempt(
                target_name,
                program_id,
                fraction,
                attempts,
            )

            corrupted = np.sort(
                rng.choice(
                    m,
                    size=k,
                    replace=False,
                )
            )

            problem = prepare_problem(
                slot_genes=slot_genes,
                corrupted_slot_idx=corrupted,
                full_program_genes=full_program_genes,
                features=features,
                target_gene_index=target["gene_index"],
            )
            mapping_rows, distances, success, _ = generate_mapping(
                problem,
                rng,
            )
            if not success:
                continue

            q = quality(distances)
            if not q["panel_valid"]:
                continue

            candidate_rows = problem["candidate"].iloc[mapping_rows]
            replacement_target_idx = candidate_rows[
                "target_corr_index"
            ].to_numpy(dtype=np.int64)

            mapped = correct_target_idx.copy()
            mapped[corrupted] = replacement_target_idx

            if len(np.unique(mapped)) != len(mapped):
                raise RuntimeError(
                    f"{target_name} {program_id} fraction={fraction}: duplicate target assignment."
                )

            edges = target_corr[
                mapped[ii],
                mapped[jj],
            ]
            rho = spearman_against_fixed(
                source_rank,
                source_norm,
                edges,
            )

            nominal_p = np.nan
            if fraction == 1.0:
                nominal_p = (
                    1
                    + int(
                        np.sum(
                            np.abs(frozen_null)
                            >= abs(rho)
                        )
                    )
                ) / (len(frozen_null) + 1)

            valid_rhos.append(rho)
            valid_rows.append(
                {
                    "target": target_name,
                    "program_id": program_id,
                    "fraction": fraction,
                    "replicate_id": len(valid_rhos),
                    "attempt_id": attempts,
                    "n_slots": m,
                    "n_corrupted_slots": k,
                    "realized_corruption_fraction": k / m,
                    "rho_edge": rho,
                    "panel_valid": 1,
                    "fraction_distance_0": q["fraction_distance_0"],
                    "fraction_distance_le1": q["fraction_distance_le1"],
                    "maximum_distance": q["maximum_distance"],
                    "mean_distance": q["mean_distance"],
                    "nominal_mapping_null_p": nominal_p,
                }
            )

        if len(valid_rhos) != REPLICATES:
            raise RuntimeError(
                f"{target_name} {program_id} fraction={fraction}: "
                f"only {len(valid_rhos)} valid replicates after {attempts} attempts."
            )

        replicate_rows.extend(valid_rows)
        vr = np.asarray(valid_rhos, dtype=np.float64)
        fraction_summaries.append(
            {
                "target": target_name,
                "program_id": program_id,
                "fraction": fraction,
                "n_replicates": REPLICATES,
                "n_corrupted_slots": k,
                "realized_corruption_fraction": k / m,
                "rho_median": float(np.median(vr)),
                "rho_q025": float(np.quantile(vr, 0.025)),
                "rho_q975": float(np.quantile(vr, 0.975)),
                "attempts": attempts,
                "invalid_attempts": attempts - REPLICATES,
            }
        )

        print(
            f"      fraction={fraction:>4.2f}; corrupted={k:4d}/{m:4d}; "
            f"median rho={np.median(vr):+.4f}; "
            f"attempts={attempts}; invalid={attempts-REPLICATES}"
        )

    rep_df = pd.DataFrame(replicate_rows)
    frac_df = pd.DataFrame(fraction_summaries)

    medians = frac_df["rho_median"].to_numpy(dtype=np.float64)
    fraction_values = frac_df["fraction"].to_numpy(dtype=np.float64)
    trend = float(
        np.corrcoef(
            rankdata(fraction_values, method="average"),
            rankdata(medians, method="average"),
        )[0, 1]
    )

    full = rep_df.loc[
        rep_df["fraction"] == 1.0
    ].copy()
    full_median = float(np.median(full["rho_edge"]))
    calibration_median_pass = bool(
        null_q025 <= full_median <= null_q975
    )

    false_specificity_rate = float(
        np.mean(
            full["nominal_mapping_null_p"].to_numpy(dtype=float)
            < ALPHA
        )
    )
    false_specificity_pass = bool(
        false_specificity_rate <= FALSE_SPECIFICITY_MAX
    )

    scope = bool(
        boolish(class_row["primary_assessable"])
        and float(class_row["rho_edge"]) > 0
        and boolish(class_row["mapping_specificity_supported"])
    )
    trend_pass = bool(trend <= TREND_GATE) if scope else None

    overall_pass = bool(
        (trend_pass if scope else True)
        and calibration_median_pass
        and false_specificity_pass
    )

    combo = {
        "target": target_name,
        "program_id": program_id,
        "primary_classification": str(class_row["primary_classification"]),
        "n_slots": m,
        "n_specificity_edges": int(len(ii)),
        "rho_zero_fraction": rho0,
        "frozen_mapping_null_q025": null_q025,
        "frozen_mapping_null_q975": null_q975,
        "trend_gate_scope": scope,
        "spearman_fraction_vs_median_rho": trend,
        "trend_gate_threshold": TREND_GATE,
        "trend_gate_pass": trend_pass,
        "full_corruption_median_rho": full_median,
        "full_corruption_median_calibration_pass": calibration_median_pass,
        "full_corruption_false_specificity_rate": false_specificity_rate,
        "false_specificity_gate_threshold": FALSE_SPECIFICITY_MAX,
        "full_corruption_false_specificity_pass": false_specificity_pass,
        "controlled_degradation_overall_pass": overall_pass,
    }

    return rep_df, frac_df, combo


def main() -> None:
    print("=" * 144)
    print("Paper 4 / TCBB - controlled mapping degradation operating-characteristic analysis")
    print("=" * 144)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Frozen scientific contract:")
    print("  Target outcomes/treatment loaded:                    NO")
    print(f"  Fractions:                                            {FRACTIONS}")
    print(f"  Replicates/nonzero fraction:                         {REPLICATES}")
    print("  Corrupted slots:                                     uniform without replacement")
    print("  Replacement genes:                                   unrelated matched background via frozen 03d generator")
    print("  Uncorrupted slots:                                   correct mapping retained")
    print("  Duplicate target assignments:                        PROHIBITED")
    print("  Primary statistic:                                   corrected 05b v2 fixed edge-subset rho")
    print("  Loading degradation endpoint:                        NO")
    print(f"  Trend gate:                                          Spearman(fraction, median rho) <= {TREND_GATE}")
    print("  100% median calibration:                             inside corrected frozen mapping-null 95% interval")
    print(f"  100% nominal false-specificity gate:                 <= {FALSE_SPECIFICITY_MAX:.0%}")
    print("=" * 144)

    for p in [
        EXEC_CONTRACT,
        OPERATING_CONTRACT,
        FINAL_CLASSIFICATION,
        SPECIFICITY_SUMMARY,
        SCANB_FEATURES,
        METABRIC_FEATURES,
        SCANB_EVAL,
        METABRIC_EVAL,
        WEIGHTS,
        SCANB_EXPR,
        SCANB_PRIMARY,
        METABRIC_EXPR,
    ]:
        require(p)

    excon = json.loads(EXEC_CONTRACT.read_text(encoding="utf-8"))
    if excon.get("status") != "FROZEN_BEFORE_CONTROLLED_DEGRADATION_RESULTS":
        raise RuntimeError("05d1 execution contract has unexpected status.")

    op = json.loads(OPERATING_CONTRACT.read_text(encoding="utf-8"))
    if op.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04a operating contract has unexpected status.")

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable.")
    device = torch.device("cuda:0")
    prop = torch.cuda.get_device_properties(0)
    print(
        f"CUDA device: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB; "
        f"capability={prop.major}.{prop.minor}"
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    classification = pd.read_csv(
        FINAL_CLASSIFICATION,
        sep="\t",
        low_memory=False,
    )
    _ = pd.read_csv(
        SPECIFICITY_SUMMARY,
        sep="\t",
        low_memory=False,
    )

    weights = pd.read_csv(
        WEIGHTS,
        sep="\t",
        dtype=str,
    ).fillna("")
    weights["Hugo_Symbol"] = weights["Hugo_Symbol"].map(canon_symbol)

    scanb_eval_df = pd.read_csv(
        SCANB_EVAL,
        sep="\t",
        dtype=str,
    ).fillna("")
    met_eval_df = pd.read_csv(
        METABRIC_EVAL,
        sep="\t",
        dtype=str,
    ).fillna("")

    scanb_eval = {
        canon_symbol(x)
        for x in scanb_eval_df.loc[
            scanb_eval_df["statistically_evaluable_corrected"].astype(str) == "1",
            "Hugo_Symbol",
        ]
    }
    met_eval = {
        canon_symbol(x)
        for x in met_eval_df.loc[
            met_eval_df["statistically_evaluable_corrected"].astype(str) == "1",
            "Hugo_Symbol",
        ]
    }

    if len(scanb_eval) != 9220:
        raise RuntimeError(f"Expected 9,220 corrected SCAN-B genes, got {len(scanb_eval)}.")
    if len(met_eval) != 8485:
        raise RuntimeError(f"Expected 8,485 corrected METABRIC genes, got {len(met_eval)}.")

    scanb_features = load_matching_features(
        "SCANB_GSE96058",
        scanb_eval,
    )
    met_features = load_matching_features(
        "METABRIC",
        met_eval,
    )

    primary = pd.read_csv(
        SCANB_PRIMARY,
        sep="\t",
        dtype=str,
    ).fillna("")
    primary_titles = list(primary["primary_title"])

    all_reps = []
    all_frac = []
    all_combo = []

    for target_name in ["SCANB_GSE96058", "METABRIC"]:
        print(f"\n[1/2] Loading corrected {target_name} target matrix ...")

        eval_set = scanb_eval if target_name == "SCANB_GSE96058" else met_eval
        target = (
            load_scanb_target(eval_set, primary_titles)
            if target_name == "SCANB_GSE96058"
            else load_metabric_target(eval_set)
        )
        features = scanb_features if target_name == "SCANB_GSE96058" else met_features

        print(
            f"  {target['z'].shape[0]:,} samples x "
            f"{target['z'].shape[1]:,} corrected-evaluable genes"
        )

        target_corr = full_corr_gpu(
            target["z"],
            device,
        )

        sub = classification.loc[
            classification["target"] == target_name
        ].copy()

        for i, row in enumerate(
            sub.to_dict(orient="records"),
            start=1,
        ):
            program_id = str(row["program_id"])

            if not boolish(row["primary_assessable"]):
                print(
                    f"\n  [{i:02d}/12] {program_id}: NOT ASSESSABLE — degradation not estimated"
                )
                continue

            print(
                f"\n  [{i:02d}/12] {program_id}: "
                f"CLASS={row['primary_classification']}; "
                f"rho0={float(row['mapping_specificity_rho']):+.4f}"
            )

            rep_df, frac_df, combo = run_combo(
                target=target,
                target_corr=target_corr,
                features=features,
                weights=weights,
                class_row=row,
            )

            all_reps.append(rep_df)
            all_frac.append(frac_df)
            all_combo.append(combo)

            print(
                f"      trend rho={combo['spearman_fraction_vs_median_rho']:+.4f} "
                f"({'PASS' if combo['trend_gate_pass'] else 'FAIL'}); "
                f"100% median={combo['full_corruption_median_rho']:+.4f} "
                f"({'PASS' if combo['full_corruption_median_calibration_pass'] else 'FAIL'}); "
                f"false-specificity={combo['full_corruption_false_specificity_rate']:.1%} "
                f"({'PASS' if combo['full_corruption_false_specificity_pass'] else 'FAIL'})"
            )

        del target_corr, target
        torch.cuda.empty_cache()
        gc.collect()

    print("\n[2/2] Writing controlled-degradation outputs ...")

    reps = pd.concat(all_reps, ignore_index=True)
    frac = pd.concat(all_frac, ignore_index=True)
    combo = pd.DataFrame(all_combo)

    reps_path = OUT_DIR / "controlled_degradation_replicates_v1.tsv"
    frac_path = OUT_DIR / "controlled_degradation_fraction_summary_v1.tsv"
    combo_path = OUT_DIR / "controlled_degradation_gate_summary_v1.tsv"

    reps.to_csv(reps_path, sep="\t", index=False)
    frac.to_csv(frac_path, sep="\t", index=False)
    combo.to_csv(combo_path, sep="\t", index=False)

    n_gate = int(combo["trend_gate_scope"].astype(bool).sum())
    n_trend_pass = int(
        combo.loc[
            combo["trend_gate_scope"].astype(bool),
            "trend_gate_pass",
        ].astype(bool).sum()
    )
    n_cal_pass = int(
        combo["full_corruption_median_calibration_pass"].astype(bool).sum()
    )
    n_false_pass = int(
        combo["full_corruption_false_specificity_pass"].astype(bool).sum()
    )
    n_overall_pass = int(
        combo["controlled_degradation_overall_pass"].astype(bool).sum()
    )

    overall_ok = (
        n_trend_pass == n_gate
        and n_cal_pass == len(combo)
        and n_false_pass == len(combo)
        and n_overall_pass == len(combo)
    )

    master = {
        "script_version": SCRIPT_VERSION,
        "status": (
            "CONTROLLED_DEGRADATION_ALL_FROZEN_GATES_PASS"
            if overall_ok
            else "CONTROLLED_DEGRADATION_ONE_OR_MORE_FROZEN_GATES_FAIL"
        ),
        "assessable_target_program_combinations": int(len(combo)),
        "formal_trend_gate_combinations": n_gate,
        "trend_gate_passes": n_trend_pass,
        "full_corruption_median_calibration_passes": n_cal_pass,
        "full_corruption_false_specificity_passes": n_false_pass,
        "overall_combination_passes": n_overall_pass,
        "fractions": FRACTIONS,
        "replicates_per_nonzero_fraction": REPLICATES,
        "primary_classification_changed": False,
        "target_outcomes_or_treatment_loaded": False,
        "replicate_file": str(reps_path),
        "fraction_summary_file": str(frac_path),
        "gate_summary_file": str(combo_path),
    }

    master_path = OUT_DIR / "controlled_mapping_degradation_v1.json"
    master_path.write_text(
        json.dumps(master, indent=2),
        encoding="utf-8",
    )

    print("\n" + "=" * 144)
    print("05d CONTROLLED MAPPING DEGRADATION: COMPLETE")
    print("=" * 144)
    print(f"Formal trend gate:                 {n_trend_pass}/{n_gate} PASS")
    print(f"100% median calibration:           {n_cal_pass}/{len(combo)} PASS")
    print(f"100% false-specificity gate:       {n_false_pass}/{len(combo)} PASS")
    print(f"Overall target/program gate:       {n_overall_pass}/{len(combo)} PASS")
    print("OVERALL STATUS:                     " + ("PASS" if overall_ok else "FAIL"))
    print()
    print(f"Replicates:       {reps_path}")
    print(f"Fraction summary: {frac_path}")
    print(f"Gate summary:     {combo_path}")
    print(f"Master JSON:      {master_path}")
    print("=" * 144)


if __name__ == "__main__":
    main()
