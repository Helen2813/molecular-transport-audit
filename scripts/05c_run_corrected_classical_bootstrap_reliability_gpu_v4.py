from __future__ import annotations

import gc
import json
import math
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import torch
except Exception as exc:
    raise RuntimeError("PyTorch is required for 05c v4.") from exc

try:
    from scipy.stats import rankdata
except Exception as exc:
    raise RuntimeError("scipy is required for exact Spearman validation.") from exc


SCRIPT_VERSION = "05c-run-corrected-classical-bootstrap-reliability-gpu-v4-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

RESTART_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_corrected_classical_bootstrap_restart_v1"
    / "corrected_classical_bootstrap_restart_v1.json"
)
EXACT_CORRECTION = (
    DATA_ROOT
    / "paper4_tcbb_exact_variation_evaluability_correction_v1"
    / "exact_variation_evaluability_correction_v1.json"
)
FINAL_CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
    / "primary_pooled_final_classification_v2.tsv"
)
DIRECT_ROOT = DATA_ROOT / "paper4_tcbb_primary_pooled_direct_preservation_v2"

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

OUT_DIR = DATA_ROOT / "paper4_tcbb_corrected_classical_bootstrap_reliability_v4"

BOOTSTRAP_DRAWS = 1000
MAX_BOOTSTRAP_ATTEMPTS = 10000
CI = 0.95
SPLIT_HALF_REPEATS = 2000

BASE_SEED = 20260916
CHECKPOINT_EVERY = 25
SPLIT_BATCH = 200

PC1_RESIDUAL_TOL = 1e-10
MAX_POWER_ITER = 400
OBSERVED_REPLAY_TOL = 5e-10
ENGINE_EQUIVALENCE_TOL = 5e-10

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


def combo_seed(target: str, program_id: str, offset: int) -> int:
    return (
        BASE_SEED
        + TARGET_INDEX[target] * 100_000
        + module_number(program_id) * 1_000
        + offset
    )


def boolish(x: object) -> bool:
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in {"1", "true", "yes", "y"}


def standardize_samples_x_genes_np(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    if not np.isfinite(x).all():
        raise RuntimeError("Non-finite values reached standardization.")
    mu = x.mean(axis=0, keepdims=True)
    sd = x.std(axis=0, ddof=1, keepdims=True)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise RuntimeError("Zero/non-finite gene SD reached standardization.")
    return (x - mu) / sd


def centered_rank_cpu(x: np.ndarray) -> tuple[np.ndarray, float]:
    r = rankdata(np.asarray(x, dtype=np.float64), method="average")
    r -= r.mean()
    ss = float(np.dot(r, r))
    if ss <= 0:
        raise RuntimeError("Degenerate rank vector.")
    return r, math.sqrt(ss)


def spearman_against_fixed_cpu(
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


def gpu_centered_average_ranks(
    values: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Exact average ranks for exact ties on CUDA.
    Ranks are 1-based, then centered. Spearman uses these centered ranks.
    """
    sorted_vals, order = torch.sort(values)
    n = sorted_vals.numel()

    new_group = torch.empty(n, dtype=torch.bool, device=values.device)
    new_group[0] = True
    if n > 1:
        new_group[1:] = sorted_vals[1:] != sorted_vals[:-1]

    group_id = torch.cumsum(new_group.to(torch.int64), dim=0) - 1
    counts = torch.bincount(group_id)
    starts0 = torch.cumsum(counts, dim=0) - counts

    avg_rank = (
        starts0.to(torch.float64)
        + (counts.to(torch.float64) - 1.0) / 2.0
        + 1.0
    )

    ranks_sorted = avg_rank[group_id]
    ranks = torch.empty(n, dtype=torch.float64, device=values.device)
    ranks[order] = ranks_sorted

    centered = ranks - (n + 1.0) / 2.0
    norm = torch.linalg.vector_norm(centered)

    del sorted_vals, order, new_group, group_id, counts, starts0, avg_rank, ranks_sorted, ranks
    return centered, norm


def gpu_spearman_with_fixed_rank(
    fixed_centered_rank_t: torch.Tensor,
    fixed_norm_t: torch.Tensor,
    values_t: torch.Tensor,
) -> torch.Tensor:
    r_t, norm_t = gpu_centered_average_ranks(values_t)
    rho = torch.dot(fixed_centered_rank_t, r_t) / (
        fixed_norm_t * norm_t
    )
    del r_t, norm_t
    return rho


def load_tcga_source(universe: pd.DataFrame) -> dict:
    raw = pd.read_csv(TCGA_EXPR, sep="\t", low_memory=False)
    row_idx = universe["source_row_index_0based"].astype(int).to_numpy()
    selected = raw.iloc[row_idx, :]

    expected = [canon_symbol(x) for x in universe["Hugo_Symbol"]]
    observed = [canon_symbol(x) for x in selected["Hugo_Symbol"]]
    if expected != observed:
        raise RuntimeError("TCGA frozen source-universe replay mismatch.")

    values = selected.iloc[:, 2:].to_numpy(dtype=np.float64)
    if np.any(values < 0) or not np.isfinite(values).all():
        raise RuntimeError("TCGA source contains invalid values.")

    z = standardize_samples_x_genes_np(np.log2(values + 1.0).T)

    return {
        "genes": expected,
        "gene_index": {g: i for i, g in enumerate(expected)},
        "z": z,
    }


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
        raise RuntimeError("No SCAN-B frozen-universe genes found.")

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
    x = x[:, eligible]

    if len(genes) != 9220:
        raise RuntimeError(
            f"Expected 9,220 corrected-evaluable SCAN-B genes, found {len(genes)}."
        )

    return {
        "name": "SCANB_GSE96058",
        "genes": genes,
        "gene_index": {g: i for i, g in enumerate(genes)},
        "x": x,
        "z": standardize_samples_x_genes_np(x),
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
        raise RuntimeError("No METABRIC frozen-universe genes found.")

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
    x = x[:, eligible]

    if len(genes) != 8485:
        raise RuntimeError(
            f"Expected 8,485 corrected-evaluable METABRIC genes, found {len(genes)}."
        )

    return {
        "name": "METABRIC",
        "genes": genes,
        "gene_index": {g: i for i, g in enumerate(genes)},
        "x": x,
        "z": standardize_samples_x_genes_np(x),
    }


def source_edge_rank_for_program(
    source: dict,
    genes: list[str],
) -> tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    cols = np.array(
        [source["gene_index"][g] for g in genes],
        dtype=np.int64,
    )
    z = source["z"][:, cols]
    corr = (z.T @ z) / (z.shape[0] - 1)
    ii, jj = np.triu_indices(len(genes), k=1)
    edges = corr[ii, jj]
    rank, norm = centered_rank_cpu(edges)
    return rank, norm, ii, jj


def count_weighted_corr(
    x_t: torch.Tensor,
    counts_t: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Exactly equivalent to an explicit ordinary bootstrap resample represented
    by sample multiplicities. The common n/(n-1) factor cancels in correlation.
    """
    total = torch.sum(counts_t)
    weights = counts_t / total

    mu = torch.sum(weights[:, None] * x_t, dim=0, keepdim=True)
    xc = x_t - mu
    var = torch.sum(weights[:, None] * xc * xc, dim=0)

    bad = (~torch.isfinite(var)) | (var <= 0)
    if bool(torch.any(bad)):
        return torch.empty(0, device=x_t.device), torch.empty(0, device=x_t.device), bad

    z = xc / torch.sqrt(var)[None, :]
    corr = z.T @ (weights[:, None] * z)
    return corr, z, bad


def leading_pc1(
    corr_t: torch.Tensor,
    z_t: torch.Tensor,
    counts_t: torch.Tensor,
    init_loading_t: torch.Tensor,
    sign_t: torch.Tensor,
    random_start_t: torch.Tensor,
) -> tuple[torch.Tensor, float, bool]:
    """
    Multi-start power iteration mirroring the 05a logic closely:
    source loading, frozen sign vector, and deterministic random start.
    """
    eps = torch.finfo(torch.float64).eps
    starts = [init_loading_t, sign_t, random_start_t]
    candidates = []

    for start_idx, start in enumerate(starts):
        v = start / torch.linalg.vector_norm(start)
        residual_value = float("inf")

        for it in range(1, MAX_POWER_ITER + 1):
            w = corr_t @ v
            wn = torch.linalg.vector_norm(w)
            if not torch.isfinite(wn) or wn <= 0:
                break
            nv = w / wn
            if torch.dot(nv, v) < 0:
                nv = -nv
            v = nv

            if it % 5 == 0 or it == MAX_POWER_ITER:
                av = corr_t @ v
                lam = torch.dot(v, av)
                residual = torch.linalg.vector_norm(
                    av - lam * v
                ) / (torch.abs(lam) + eps)
                residual_value = float(residual.detach().cpu())
                if residual_value <= PC1_RESIDUAL_TOL:
                    break

        av = corr_t @ v
        lam = float(torch.dot(v, av).detach().cpu())
        candidates.append(
            (
                lam,
                residual_value,
                it,
                start_idx,
                v.detach().clone(),
            )
        )

    candidates.sort(key=lambda x: x[0], reverse=True)
    lam, residual, iters, chosen_start, v = candidates[0]

    fallback = False
    if residual > PC1_RESIDUAL_TOL:
        evals, evecs = torch.linalg.eigh(corr_t)
        v = evecs[:, -1]
        lam_t = evals[-1]
        av = corr_t @ v
        residual_t = torch.linalg.vector_norm(
            av - lam_t * v
        ) / (torch.abs(lam_t) + eps)
        residual = float(residual_t.detach().cpu())
        fallback = True

    weights = counts_t / torch.sum(counts_t)
    scores = z_t @ v
    sign_score = torch.mean(z_t * sign_t[None, :], dim=1)

    score_mean = torch.sum(weights * scores)
    sign_mean = torch.sum(weights * sign_score)
    sc = scores - score_mean
    sg = sign_score - sign_mean

    num = torch.sum(weights * sc * sg)
    den = torch.sqrt(
        torch.sum(weights * sc * sc)
        * torch.sum(weights * sg * sg)
    )
    if den <= 0 or not torch.isfinite(den):
        raise RuntimeError("PC1 orientation statistic is degenerate.")

    if num / den < 0:
        v = -v

    return v, residual, fallback


def explicit_resample_corr_reference(
    x: np.ndarray,
    idx: np.ndarray,
    k: int = 80,
) -> np.ndarray:
    xb = x[idx, :k]
    z = standardize_samples_x_genes_np(xb)
    return (z.T @ z) / (z.shape[0] - 1)


def reliability_diagnostics(
    target_z: np.ndarray,
    signs: np.ndarray,
    repeats: int,
    seed: int,
    device: torch.device,
) -> dict:
    n, m = target_z.shape
    zt = torch.as_tensor(
        target_z,
        dtype=torch.float64,
        device=device,
    )
    st = torch.as_tensor(
        signs,
        dtype=torch.float64,
        device=device,
    )
    signed = zt * st[None, :]

    n1 = m // 2
    n2 = m - n1

    rng = np.random.default_rng(seed)
    raw = np.empty(repeats, dtype=np.float64)
    sb = np.empty(repeats, dtype=np.float64)

    done = 0
    while done < repeats:
        k = min(SPLIT_BATCH, repeats - done)
        mask = np.zeros((m, k), dtype=np.float64)

        for col in range(k):
            idx = rng.permutation(m)[:n1]
            mask[idx, col] = 1.0

        mt = torch.as_tensor(mask, dtype=torch.float64, device=device)
        sum1 = signed @ mt
        total = torch.sum(signed, dim=1, keepdim=True)
        sum2 = total - sum1

        s1 = sum1 / n1
        s2 = sum2 / n2
        s1 = s1 - s1.mean(dim=0, keepdim=True)
        s2 = s2 - s2.mean(dim=0, keepdim=True)

        num = torch.sum(s1 * s2, dim=0)
        den = torch.sqrt(
            torch.sum(s1 * s1, dim=0)
            * torch.sum(s2 * s2, dim=0)
        )
        r = (num / den).detach().cpu().numpy()

        raw[done : done + k] = r
        sb[done : done + k] = 2.0 * r / (1.0 + r)

        del mt, sum1, total, sum2, s1, s2, num, den
        done += k

    full = torch.mean(signed, dim=1)
    loo = (m * full[:, None] - signed) / (m - 1)

    full_c = full - full.mean()
    loo_c = loo - loo.mean(dim=0, keepdim=True)
    num = torch.sum(full_c[:, None] * loo_c, dim=0)
    den = torch.linalg.vector_norm(full_c) * torch.sqrt(
        torch.sum(loo_c * loo_c, dim=0)
    )
    loo_cor = (num / den).detach().cpu().numpy()

    del zt, st, signed, full, loo, full_c, loo_c, num, den
    torch.cuda.empty_cache()

    return {
        "split_half_raw": raw,
        "split_half_spearman_brown": sb,
        "leave_one_out_correlations": loo_cor,
    }


def output_paths(target: str, program_id: str) -> dict[str, Path]:
    d = OUT_DIR / "per_program" / target
    d.mkdir(parents=True, exist_ok=True)
    return {
        "checkpoint": d / f"{program_id}_classical_bootstrap_checkpoint_v4.npz",
        "arrays": d / f"{program_id}_classical_bootstrap_reliability_arrays_v4.npz",
        "result": d / f"{program_id}_classical_bootstrap_reliability_result_v4.json",
        "degenerate": d / f"{program_id}_degenerate_bootstrap_gene_counts_v4.tsv",
    }


def run_program(
    target: dict,
    source: dict,
    row: dict,
    weights_df: pd.DataFrame,
    device: torch.device,
) -> dict:
    target_name = str(row["target"])
    program_id = str(row["program_id"])
    paths = output_paths(target_name, program_id)

    if paths["result"].exists() and paths["arrays"].exists():
        try:
            old = json.loads(paths["result"].read_text(encoding="utf-8"))
        except Exception:
            old = {}
        if (
            old.get("script_version") == SCRIPT_VERSION
            and int(old.get("valid_bootstrap_draws", -1)) == BOOTSTRAP_DRAWS
            and int(old.get("split_half_repeats", -1)) == SPLIT_HALF_REPEATS
        ):
            print("    deterministic checkpoint found — reusing completed result")
            return old

    gene_file = (
        DIRECT_ROOT
        / "per_program"
        / target_name
        / f"{program_id}_evaluable_genes_v2.tsv"
    )
    require(gene_file)

    eg = pd.read_csv(gene_file, sep="\t", dtype=str).fillna("")
    genes = [canon_symbol(x) for x in eg["Hugo_Symbol"]]
    m = len(genes)

    tcols = np.array(
        [target["gene_index"][g] for g in genes],
        dtype=np.int64,
    )
    x = target["x"][:, tcols]
    z_observed = target["z"][:, tcols]

    wdf = weights_df.loc[
        weights_df["program_id"] == program_id
    ].copy()
    wdf["Hugo_Symbol"] = wdf["Hugo_Symbol"].map(canon_symbol)
    wdf = wdf.set_index("Hugo_Symbol", drop=False)

    source_loading = np.array(
        [
            float(wdf.loc[g, "source_pc1_loading"])
            for g in genes
        ],
        dtype=np.float64,
    )
    source_sign = np.array(
        [
            float(wdf.loc[g, "source_pc1_loading_sign"])
            for g in genes
        ],
        dtype=np.float64,
    )

    source_loading_rank, source_loading_norm = centered_rank_cpu(
        source_loading
    )
    source_edge_rank, source_edge_norm, tri_i, tri_j = (
        source_edge_rank_for_program(source, genes)
    )

    # GPU fixed objects.
    x_t = torch.as_tensor(x, dtype=torch.float64, device=device)
    source_edge_rank_t = torch.as_tensor(
        source_edge_rank,
        dtype=torch.float64,
        device=device,
    )
    source_edge_norm_t = torch.tensor(
        source_edge_norm,
        dtype=torch.float64,
        device=device,
    )
    source_loading_rank_t = torch.as_tensor(
        source_loading_rank,
        dtype=torch.float64,
        device=device,
    )
    source_loading_norm_t = torch.tensor(
        source_loading_norm,
        dtype=torch.float64,
        device=device,
    )
    init_t = torch.as_tensor(
        source_loading,
        dtype=torch.float64,
        device=device,
    )
    sign_t = torch.as_tensor(
        source_sign,
        dtype=torch.float64,
        device=device,
    )
    ii_t = torch.as_tensor(tri_i, dtype=torch.long, device=device)
    jj_t = torch.as_tensor(tri_j, dtype=torch.long, device=device)

    rng_start = np.random.default_rng(
        combo_seed(target_name, program_id, 311)
    )
    random_start_t = torch.as_tensor(
        rng_start.normal(size=m),
        dtype=torch.float64,
        device=device,
    )

    # -----------------------------------------------------------------
    # Mandatory observed replay using count=1 for every original sample.
    # -----------------------------------------------------------------
    ones_t = torch.ones(
        x.shape[0],
        dtype=torch.float64,
        device=device,
    )
    corr_obs, z_obs_t, bad_obs = count_weighted_corr(x_t, ones_t)
    if corr_obs.numel() == 0 or bool(torch.any(bad_obs)):
        raise RuntimeError(
            f"{target_name} {program_id}: corrected observed gene set is still degenerate."
        )

    target_edges_obs = corr_obs[ii_t, jj_t]
    replay_edge_t = gpu_spearman_with_fixed_rank(
        source_edge_rank_t,
        source_edge_norm_t,
        target_edges_obs,
    )
    replay_edge = float(replay_edge_t.detach().cpu())

    load_obs_t, replay_resid, replay_fallback = leading_pc1(
        corr_obs,
        z_obs_t,
        ones_t,
        init_t,
        sign_t,
        random_start_t,
    )
    replay_load_t = gpu_spearman_with_fixed_rank(
        source_loading_rank_t,
        source_loading_norm_t,
        load_obs_t,
    )
    replay_load = float(replay_load_t.detach().cpu())

    observed_edge = float(row["rho_edge"])
    observed_load = float(row["rho_load"])

    if abs(replay_edge - observed_edge) > OBSERVED_REPLAY_TOL:
        raise RuntimeError(
            f"{target_name} {program_id}: observed edge replay mismatch "
            f"{replay_edge} vs {observed_edge}"
        )
    if abs(replay_load - observed_load) > OBSERVED_REPLAY_TOL:
        raise RuntimeError(
            f"{target_name} {program_id}: observed loading replay mismatch "
            f"{replay_load} vs {observed_load}"
        )

    # GPU-rank vs scipy-rank replay check on observed target edges.
    target_edges_obs_np = target_edges_obs.detach().cpu().numpy()
    replay_edge_cpu = spearman_against_fixed_cpu(
        source_edge_rank,
        source_edge_norm,
        target_edges_obs_np,
    )
    if abs(replay_edge_cpu - replay_edge) > 2e-12:
        raise RuntimeError(
            f"{target_name} {program_id}: GPU average-rank Spearman "
            f"differs from scipy: {replay_edge} vs {replay_edge_cpu}"
        )

    del (
        corr_obs,
        z_obs_t,
        bad_obs,
        target_edges_obs,
        replay_edge_t,
        load_obs_t,
        replay_load_t,
        target_edges_obs_np,
        ones_t,
    )
    torch.cuda.empty_cache()

    # Reliability is descriptive and uses the corrected full pooled target.
    rel = reliability_diagnostics(
        target_z=z_observed,
        signs=source_sign,
        repeats=SPLIT_HALF_REPEATS,
        seed=combo_seed(target_name, program_id, 401),
        device=device,
    )

    # -----------------------------------------------------------------
    # Bootstrap state. Arrays index VALID draws.
    # -----------------------------------------------------------------
    edge_rho = np.full(BOOTSTRAP_DRAWS, np.nan, dtype=np.float64)
    load_rho = np.full(BOOTSTRAP_DRAWS, np.nan, dtype=np.float64)
    pc1_resid = np.full(BOOTSTRAP_DRAWS, np.nan, dtype=np.float64)
    pc1_fallback = np.zeros(BOOTSTRAP_DRAWS, dtype=np.int8)

    completed = 0
    attempts = 0
    invalid_attempts = 0
    degenerate_gene_counts = np.zeros(m, dtype=np.int64)

    seed = combo_seed(target_name, program_id, 301)

    if paths["checkpoint"].exists():
        ck = np.load(paths["checkpoint"])
        if (
            int(ck["seed"]) == seed
            and int(ck["draws"]) == BOOTSTRAP_DRAWS
            and int(ck["max_attempts"]) == MAX_BOOTSTRAP_ATTEMPTS
        ):
            completed = int(ck["completed"])
            attempts = int(ck["attempts"])
            invalid_attempts = int(ck["invalid_attempts"])
            edge_rho[:] = ck["edge_rho"]
            load_rho[:] = ck["load_rho"]
            pc1_resid[:] = ck["pc1_resid"]
            pc1_fallback[:] = ck["pc1_fallback"]
            degenerate_gene_counts[:] = ck["degenerate_gene_counts"]
            print(
                f"    resuming valid={completed}/{BOOTSTRAP_DRAWS}; "
                f"attempts={attempts}; invalid={invalid_attempts}"
            )

    rng = np.random.default_rng(seed)

    # Replay RNG by ATTEMPT count.
    for _ in range(attempts):
        _ = rng.integers(
            0,
            x.shape[0],
            size=x.shape[0],
            endpoint=False,
        )

    engine_validated = False
    t0 = time.perf_counter()
    last_checkpoint = completed

    while (
        completed < BOOTSTRAP_DRAWS
        and attempts < MAX_BOOTSTRAP_ATTEMPTS
    ):
        idx = rng.integers(
            0,
            x.shape[0],
            size=x.shape[0],
            endpoint=False,
        )
        attempts += 1

        counts_np = np.bincount(
            idx,
            minlength=x.shape[0],
        ).astype(np.float64)
        counts_t = torch.as_tensor(
            counts_np,
            dtype=torch.float64,
            device=device,
        )

        corr, z_t, bad_t = count_weighted_corr(
            x_t,
            counts_t,
        )

        if corr.numel() == 0:
            bad = bad_t.detach().cpu().numpy().astype(bool)
            degenerate_gene_counts[bad] += 1
            invalid_attempts += 1
            del counts_t, corr, z_t, bad_t
            continue

        # First valid draw: prove count-weight implementation equals an
        # explicit sample-with-replacement calculation on a fixed gene subset.
        if not engine_validated:
            k = min(80, m)
            ref = explicit_resample_corr_reference(
                x,
                idx,
                k=k,
            )
            got = corr[:k, :k].detach().cpu().numpy()
            diff = float(np.max(np.abs(ref - got)))
            if diff > ENGINE_EQUIVALENCE_TOL:
                raise RuntimeError(
                    f"{target_name} {program_id}: count-weight bootstrap engine "
                    f"does not match explicit resampling; max|Δ|={diff:.3e}"
                )
            print(
                f"    count-weight ↔ explicit-resample validation "
                f"max|Δ|={diff:.3e}"
            )
            engine_validated = True

        target_edges_t = corr[ii_t, jj_t]
        erho_t = gpu_spearman_with_fixed_rank(
            source_edge_rank_t,
            source_edge_norm_t,
            target_edges_t,
        )
        edge_rho[completed] = float(erho_t.detach().cpu())

        load_t, residual, fallback = leading_pc1(
            corr,
            z_t,
            counts_t,
            init_t,
            sign_t,
            random_start_t,
        )
        lrho_t = gpu_spearman_with_fixed_rank(
            source_loading_rank_t,
            source_loading_norm_t,
            load_t,
        )
        load_rho[completed] = float(lrho_t.detach().cpu())
        pc1_resid[completed] = residual
        pc1_fallback[completed] = int(fallback)

        completed += 1

        del (
            counts_t,
            corr,
            z_t,
            bad_t,
            target_edges_t,
            erho_t,
            load_t,
            lrho_t,
        )

        if (
            completed - last_checkpoint >= CHECKPOINT_EVERY
            or completed == BOOTSTRAP_DRAWS
        ):
            np.savez_compressed(
                paths["checkpoint"],
                seed=np.array(seed, dtype=np.int64),
                draws=np.array(BOOTSTRAP_DRAWS, dtype=np.int64),
                max_attempts=np.array(MAX_BOOTSTRAP_ATTEMPTS, dtype=np.int64),
                completed=np.array(completed, dtype=np.int64),
                attempts=np.array(attempts, dtype=np.int64),
                invalid_attempts=np.array(invalid_attempts, dtype=np.int64),
                edge_rho=edge_rho,
                load_rho=load_rho,
                pc1_resid=pc1_resid,
                pc1_fallback=pc1_fallback,
                degenerate_gene_counts=degenerate_gene_counts,
            )
            last_checkpoint = completed

            print(
                f"    valid {completed:4d}/{BOOTSTRAP_DRAWS}; "
                f"attempts={attempts:4d}; invalid={invalid_attempts:3d}; "
                f"edge median={np.nanmedian(edge_rho[:completed]):+.4f}; "
                f"load median={np.nanmedian(load_rho[:completed]):+.4f}; "
                f"elapsed={time.perf_counter()-t0:.1f}s"
            )

    if completed != BOOTSTRAP_DRAWS:
        raise RuntimeError(
            f"{target_name} {program_id}: only {completed} valid draws in "
            f"{attempts} attempts; frozen cap={MAX_BOOTSTRAP_ATTEMPTS}."
        )

    deg_df = pd.DataFrame(
        {
            "Hugo_Symbol": genes,
            "degenerate_attempt_count": degenerate_gene_counts,
        }
    ).sort_values(
        "degenerate_attempt_count",
        ascending=False,
    )
    deg_df.to_csv(
        paths["degenerate"],
        sep="\t",
        index=False,
    )

    alpha = (1.0 - CI) / 2.0

    result = {
        "script_version": SCRIPT_VERSION,
        "target": target_name,
        "program_id": program_id,
        "n_samples": int(x.shape[0]),
        "n_corrected_evaluable_genes": m,
        "n_edges": int(len(tri_i)),
        "valid_bootstrap_draws": BOOTSTRAP_DRAWS,
        "bootstrap_attempts": attempts,
        "invalid_bootstrap_attempts": invalid_attempts,
        "valid_attempt_fraction": float(
            BOOTSTRAP_DRAWS / attempts
        ),
        "bootstrap_seed": seed,
        "observed_rho_edge": observed_edge,
        "observed_rho_load": observed_load,
        "observed_replay_rho_edge": replay_edge,
        "observed_replay_rho_load": replay_load,
        "observed_gpu_vs_scipy_edge_spearman_absdiff": abs(
            replay_edge - replay_edge_cpu
        ),
        "observed_pc1_residual": replay_resid,
        "observed_pc1_exact_eigh_fallback": int(replay_fallback),
        "edge_bootstrap_median": float(np.median(edge_rho)),
        "edge_bootstrap_ci_low": float(np.quantile(edge_rho, alpha)),
        "edge_bootstrap_ci_high": float(np.quantile(edge_rho, 1.0-alpha)),
        "loading_bootstrap_median": float(np.median(load_rho)),
        "loading_bootstrap_ci_low": float(np.quantile(load_rho, alpha)),
        "loading_bootstrap_ci_high": float(np.quantile(load_rho, 1.0-alpha)),
        "loading_exact_eigh_fallbacks": int(pc1_fallback.sum()),
        "loading_max_pc1_residual": float(np.max(pc1_resid)),
        "genes_ever_degenerate": int(
            np.sum(degenerate_gene_counts > 0)
        ),
        "maximum_single_gene_degenerate_attempts": int(
            np.max(degenerate_gene_counts)
        ),
        "split_half_repeats": SPLIT_HALF_REPEATS,
        "split_half_raw_median": float(np.median(rel["split_half_raw"])),
        "split_half_raw_q05": float(np.quantile(rel["split_half_raw"], 0.05)),
        "split_half_raw_q95": float(np.quantile(rel["split_half_raw"], 0.95)),
        "split_half_spearman_brown_median": float(
            np.median(rel["split_half_spearman_brown"])
        ),
        "split_half_spearman_brown_q05": float(
            np.quantile(rel["split_half_spearman_brown"], 0.05)
        ),
        "split_half_spearman_brown_q95": float(
            np.quantile(rel["split_half_spearman_brown"], 0.95)
        ),
        "leave_one_out_min": float(
            np.min(rel["leave_one_out_correlations"])
        ),
        "leave_one_out_median": float(
            np.median(rel["leave_one_out_correlations"])
        ),
        "elapsed_bootstrap_seconds": time.perf_counter() - t0,
        "primary_classification_changed": False,
    }

    np.savez_compressed(
        paths["arrays"],
        edge_bootstrap_rho=edge_rho,
        loading_bootstrap_rho=load_rho,
        loading_pc1_residual=pc1_resid,
        loading_exact_eigh_fallback=pc1_fallback,
        split_half_raw=rel["split_half_raw"],
        split_half_spearman_brown=rel["split_half_spearman_brown"],
        leave_one_out_correlations=rel["leave_one_out_correlations"],
        degenerate_gene_counts=degenerate_gene_counts,
    )

    paths["result"].write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )

    del (
        x_t,
        source_edge_rank_t,
        source_edge_norm_t,
        source_loading_rank_t,
        source_loading_norm_t,
        init_t,
        sign_t,
        ii_t,
        jj_t,
        random_start_t,
    )
    torch.cuda.empty_cache()

    return result


def main() -> None:
    print("=" * 144)
    print("Paper 4 / TCBB - corrected ORIGINAL classical target bootstrap + reliability diagnostics")
    print("=" * 144)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Frozen scientific contract:")
    print("  Target outcomes/treatment loaded:                    NO")
    print(f"  Valid ordinary bootstrap draws:                      {BOOTSTRAP_DRAWS:,}")
    print(f"  Maximum attempts:                                    {MAX_BOOTSTRAP_ATTEMPTS:,}")
    print("  Corrected gene eligibility:                          finite AND exact max>min")
    print("  Source object:                                       FIXED")
    print("  Module gene set within target/program:               FIXED")
    print("  Edge statistic:                                      full corrected-evaluable-edge rho_edge")
    print("  Loading statistic:                                   recomputed target PC1 -> rho_load")
    print(f"  Percentile interval:                                 {int(CI*100)}%")
    print("  Count-weight computation:                            exact equivalent of sample resampling")
    print(f"  Split-half repeats:                                  {SPLIT_HALF_REPEATS:,}")
    print("  Bayesian bootstrap:                                  NOT USED")
    print("  Primary classification:                              corrected 05b v2, NOT CHANGED")
    print("=" * 144)

    for p in [
        RESTART_CONTRACT,
        EXACT_CORRECTION,
        FINAL_CLASSIFICATION,
        UNIVERSE,
        WEIGHTS,
        TCGA_EXPR,
        SCANB_EXPR,
        SCANB_PRIMARY,
        METABRIC_EXPR,
    ]:
        require(p)

    restart = json.loads(
        RESTART_CONTRACT.read_text(encoding="utf-8")
    )
    if restart.get("status") != "FROZEN_BEFORE_CORRECTED_CLASSICAL_BOOTSTRAP_RESULTS":
        raise RuntimeError("05c3 restart contract has unexpected status.")

    exact = json.loads(
        EXACT_CORRECTION.read_text(encoding="utf-8")
    )
    if exact.get("status") != "FROZEN_CORRECTION_BEFORE_RECOMPUTED_PRESERVATION":
        raise RuntimeError("04i3 exact-variation correction has unexpected status.")

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

    universe = pd.read_csv(
        UNIVERSE,
        sep="\t",
        dtype=str,
    ).fillna("")
    universe["source_gene_rank_by_MAD"] = pd.to_numeric(
        universe["source_gene_rank_by_MAD"],
        errors="raise",
    ).astype(int)
    universe = universe.sort_values(
        "source_gene_rank_by_MAD"
    ).reset_index(drop=True)

    if len(universe) != 10_000:
        raise RuntimeError(
            f"Expected 10,000 source genes, found {len(universe)}."
        )

    universe_set = {
        canon_symbol(x)
        for x in universe["Hugo_Symbol"]
    }

    weights_df = pd.read_csv(
        WEIGHTS,
        sep="\t",
        dtype=str,
    ).fillna("")

    primary = pd.read_csv(
        SCANB_PRIMARY,
        sep="\t",
        dtype=str,
    ).fillna("")
    primary_titles = list(primary["primary_title"])
    if len(primary_titles) != 3273:
        raise RuntimeError(
            f"Expected 3,273 SCAN-B primary samples, found {len(primary_titles)}."
        )

    print("\n[1/4] Replaying frozen TCGA source ...")
    source = load_tcga_source(universe)
    print(
        f"  {source['z'].shape[0]:,} samples x "
        f"{source['z'].shape[1]:,} genes"
    )

    results = []

    for target_name in [
        "SCANB_GSE96058",
        "METABRIC",
    ]:
        print(f"\n[2/4] Loading corrected {target_name} ...")
        target = (
            load_scanb_target(
                universe_set,
                primary_titles,
            )
            if target_name == "SCANB_GSE96058"
            else load_metabric_target(universe_set)
        )
        print(
            f"  {target['x'].shape[0]:,} samples x "
            f"{target['x'].shape[1]:,} corrected-evaluable genes"
        )

        print("\n" + "=" * 144)
        print(
            f"[3/4] CORRECTED CLASSICAL BOOTSTRAP + RELIABILITY — {target_name}"
        )
        print("=" * 144)

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
                    f"\n  [{i:02d}/12] {program_id}: "
                    "NOT ASSESSABLE — uncertainty not estimated"
                )
                results.append(
                    {
                        "target": target_name,
                        "program_id": program_id,
                        "primary_assessable": False,
                        "reason": "corrected frozen assessability guard failed",
                    }
                )
                continue

            print(
                f"\n  [{i:02d}/12] {program_id}: "
                f"CLASS={row['primary_classification']}; "
                f"edge={float(row['rho_edge']):+.4f}; "
                f"load={float(row['rho_load']):+.4f}"
            )

            result = run_program(
                target=target,
                source=source,
                row=row,
                weights_df=weights_df,
                device=device,
            )
            results.append(result)

            print(
                f"    edge 95% CI="
                f"[{result['edge_bootstrap_ci_low']:+.4f},"
                f"{result['edge_bootstrap_ci_high']:+.4f}]; "
                f"loading CI="
                f"[{result['loading_bootstrap_ci_low']:+.4f},"
                f"{result['loading_bootstrap_ci_high']:+.4f}]"
            )
            print(
                f"    split-half raw median="
                f"{result['split_half_raw_median']:.4f}; "
                f"Spearman-Brown="
                f"{result['split_half_spearman_brown_median']:.4f}; "
                f"LOO min={result['leave_one_out_min']:.6f}"
            )

        del target
        torch.cuda.empty_cache()
        gc.collect()

    print("\n[4/4] Writing corrected classical-bootstrap/reliability summary ...")

    out_df = pd.DataFrame(results)
    summary_tsv = (
        OUT_DIR
        / "primary_pooled_corrected_classical_bootstrap_reliability_summary_v4.tsv"
    )
    out_df.to_csv(
        summary_tsv,
        sep="\t",
        index=False,
    )

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "CORRECTED_PRIMARY_POOLED_CLASSICAL_BOOTSTRAP_RELIABILITY_COMPLETE",
        "valid_bootstrap_draws": BOOTSTRAP_DRAWS,
        "maximum_attempts": MAX_BOOTSTRAP_ATTEMPTS,
        "interval": CI,
        "split_half_repeats": SPLIT_HALF_REPEATS,
        "bayesian_bootstrap_used": False,
        "primary_classification_changed": False,
        "summary_file": str(summary_tsv),
    }

    master_json = (
        OUT_DIR
        / "primary_pooled_corrected_classical_bootstrap_reliability_v4.json"
    )
    master_json.write_text(
        json.dumps(master, indent=2),
        encoding="utf-8",
    )

    print("\n" + "=" * 144)
    print("05c v4 CORRECTED CLASSICAL BOOTSTRAP / RELIABILITY: COMPLETE")
    print("=" * 144)
    print("Primary classifications were NOT changed.")
    print("Bayesian bootstrap was NOT used.")
    print("No target outcomes or treatment variables were loaded.")
    print()
    print(f"Summary TSV: {summary_tsv}")
    print(f"Master JSON: {master_json}")
    print("=" * 144)


if __name__ == "__main__":
    main()
