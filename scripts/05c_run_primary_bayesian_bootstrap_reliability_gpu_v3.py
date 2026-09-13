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
    raise RuntimeError("PyTorch is required for 05c v3.") from exc

try:
    from scipy.stats import rankdata
except Exception as exc:
    raise RuntimeError("scipy is required for exact Spearman statistics.") from exc


SCRIPT_VERSION = "05c-run-primary-bayesian-bootstrap-reliability-gpu-v3-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

BAYESIAN_AMENDMENT = (
    DATA_ROOT
    / "paper4_tcbb_bayesian_bootstrap_amendment_v1"
    / "bayesian_bootstrap_uncertainty_amendment_v1.json"
)
FINAL_CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v1"
    / "primary_pooled_final_classification_v1.tsv"
)
DIRECT_ROOT = DATA_ROOT / "paper4_tcbb_primary_pooled_direct_preservation_v1"

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

OUT_DIR = DATA_ROOT / "paper4_tcbb_primary_bayesian_bootstrap_reliability_v3"

BAYESIAN_DRAWS = 1000
INTERVAL = 0.95
SPLIT_HALF_REPEATS = 2000
BASE_SEED = 20260916
CHECKPOINT_EVERY = 25
SPLIT_BATCH = 200
PC1_RESIDUAL_TOL = 1e-10
MAX_POWER_ITER = 400
REPLAY_TOL = 5e-10

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


def centered_rank(x: np.ndarray) -> tuple[np.ndarray, float]:
    r = rankdata(np.asarray(x, dtype=np.float64), method="average")
    r -= r.mean()
    ss = float(np.dot(r, r))
    if ss <= 0:
        raise RuntimeError("Degenerate rank vector.")
    return r, math.sqrt(ss)


def spearman_against_fixed_rank(
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


def load_tcga_source(universe: pd.DataFrame) -> dict:
    raw = pd.read_csv(TCGA_EXPR, sep="\t", low_memory=False)
    row_idx = universe["source_row_index_0based"].astype(int).to_numpy()
    selected = raw.iloc[row_idx, :]

    expected = [canon_symbol(x) for x in universe["Hugo_Symbol"]]
    observed = [canon_symbol(x) for x in selected["Hugo_Symbol"]]
    if expected != observed:
        raise RuntimeError("TCGA frozen source-universe replay mismatch.")

    values = selected.iloc[:, 2:].to_numpy(dtype=np.float64)
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

    df = pd.concat(parts, ignore_index=True)
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    x = df.to_numpy(dtype=np.float64).T
    finite = np.isfinite(x).all(axis=0)
    sd = np.full(x.shape[1], np.nan)
    sd[finite] = x[:, finite].std(axis=0, ddof=1)
    eligible = finite & np.isfinite(sd) & (sd > 0)

    genes_all = list(df.index)
    genes = [g for g, ok in zip(genes_all, eligible) if ok]
    x = x[:, eligible]

    if len(genes) != 9225:
        raise RuntimeError(f"Expected 9,225 SCAN-B evaluable genes, got {len(genes)}.")

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

    df = pd.concat(parts, ignore_index=True)
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    x = df.to_numpy(dtype=np.float64).T
    finite = np.isfinite(x).all(axis=0)
    sd = np.full(x.shape[1], np.nan)
    sd[finite] = x[:, finite].std(axis=0, ddof=1)
    eligible = finite & np.isfinite(sd) & (sd > 0)

    genes_all = list(df.index)
    genes = [g for g, ok in zip(genes_all, eligible) if ok]
    x = x[:, eligible]

    if len(genes) != 8485:
        raise RuntimeError(f"Expected 8,485 METABRIC evaluable genes, got {len(genes)}.")

    return {
        "name": "METABRIC",
        "genes": genes,
        "gene_index": {g: i for i, g in enumerate(genes)},
        "x": x,
        "z": standardize_samples_x_genes_np(x),
    }


def source_edge_rank(
    source: dict,
    genes: list[str],
) -> tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    cols = np.array([source["gene_index"][g] for g in genes], dtype=np.int64)
    z = source["z"][:, cols]
    corr = (z.T @ z) / (z.shape[0] - 1)
    ii, jj = np.triu_indices(len(genes), k=1)
    edges = corr[ii, jj]
    return (*centered_rank(edges), ii, jj)


def leading_pc1_weighted(
    corr_t: torch.Tensor,
    z_t: torch.Tensor,
    weights_t: torch.Tensor,
    init_t: torch.Tensor,
    sign_t: torch.Tensor,
) -> tuple[torch.Tensor, float, bool]:
    eps = torch.finfo(torch.float64).eps
    v = init_t / torch.linalg.vector_norm(init_t)

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
            residual = torch.linalg.vector_norm(av - lam * v) / (
                torch.abs(lam) + eps
            )
            residual_value = float(residual.detach().cpu())
            if residual_value <= PC1_RESIDUAL_TOL:
                break

    fallback = False
    if residual_value > PC1_RESIDUAL_TOL:
        evals, evecs = torch.linalg.eigh(corr_t)
        v = evecs[:, -1]
        lam = evals[-1]
        av = corr_t @ v
        residual = torch.linalg.vector_norm(av - lam * v) / (
            torch.abs(lam) + eps
        )
        residual_value = float(residual.detach().cpu())
        fallback = True

    scores = z_t @ v
    sign_score = torch.mean(z_t * sign_t[None, :], dim=1)

    score_mean = torch.sum(weights_t * scores)
    sign_mean = torch.sum(weights_t * sign_score)

    sc = scores - score_mean
    sg = sign_score - sign_mean

    num = torch.sum(weights_t * sc * sg)
    den = torch.sqrt(
        torch.sum(weights_t * sc * sc)
        * torch.sum(weights_t * sg * sg)
    )
    if den <= 0:
        raise RuntimeError("Degenerate weighted PC1 orientation score.")

    if num / den < 0:
        v = -v

    return v, residual_value, fallback


def weighted_corr(
    x_t: torch.Tensor,
    weights_t: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor]:
    mu = torch.sum(weights_t[:, None] * x_t, dim=0, keepdim=True)
    xc = x_t - mu
    var = torch.sum(weights_t[:, None] * xc * xc, dim=0)

    if torch.any(~torch.isfinite(var)) or torch.any(var <= 0):
        raise RuntimeError(
            "Bayesian-bootstrap weighted variance became nonfinite/nonpositive."
        )

    z = xc / torch.sqrt(var)[None, :]
    corr = z.T @ (weights_t[:, None] * z)
    return corr, z


def reliability_diagnostics(
    target_z: np.ndarray,
    signs: np.ndarray,
    repeats: int,
    seed: int,
    device: torch.device,
) -> dict:
    n, m = target_z.shape
    zt = torch.as_tensor(target_z, dtype=torch.float64, device=device)
    st = torch.as_tensor(signs, dtype=torch.float64, device=device)
    signed = zt * st[None, :]

    n1 = m // 2
    n2 = m - n1

    rng = np.random.default_rng(seed)
    raw = np.empty(repeats, dtype=np.float64)
    sb = np.empty(repeats, dtype=np.float64)

    done = 0
    while done < repeats:
        k = min(SPLIT_BATCH, repeats - done)
        w1 = np.zeros((m, k), dtype=np.float64)

        for col in range(k):
            idx = rng.permutation(m)[:n1]
            w1[idx, col] = 1.0

        w1t = torch.as_tensor(w1, dtype=torch.float64, device=device)
        sum1 = signed @ w1t
        total = torch.sum(signed, dim=1, keepdim=True)
        sum2 = total - sum1

        score1 = sum1 / n1
        score2 = sum2 / n2

        score1 = score1 - score1.mean(dim=0, keepdim=True)
        score2 = score2 - score2.mean(dim=0, keepdim=True)

        num = torch.sum(score1 * score2, dim=0)
        den = torch.sqrt(
            torch.sum(score1 * score1, dim=0)
            * torch.sum(score2 * score2, dim=0)
        )
        r = (num / den).detach().cpu().numpy()

        raw[done : done + k] = r
        sb[done : done + k] = 2.0 * r / (1.0 + r)

        del w1t, sum1, sum2, score1, score2, num, den
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


def paths_for(target: str, program_id: str) -> dict[str, Path]:
    d = OUT_DIR / "per_program" / target
    d.mkdir(parents=True, exist_ok=True)
    return {
        "checkpoint": d / f"{program_id}_bayesian_bootstrap_checkpoint_v3.npz",
        "arrays": d / f"{program_id}_bayesian_bootstrap_reliability_arrays_v3.npz",
        "result": d / f"{program_id}_bayesian_bootstrap_reliability_result_v3.json",
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
    paths = paths_for(target_name, program_id)

    if paths["result"].exists() and paths["arrays"].exists():
        try:
            old = json.loads(paths["result"].read_text(encoding="utf-8"))
        except Exception:
            old = {}
        if (
            old.get("script_version") == SCRIPT_VERSION
            and int(old.get("bayesian_bootstrap_draws", -1)) == BAYESIAN_DRAWS
            and int(old.get("split_half_repeats", -1)) == SPLIT_HALF_REPEATS
        ):
            print("    deterministic checkpoint found — reusing completed result")
            return old

    gene_file = (
        DIRECT_ROOT
        / "per_program"
        / target_name
        / f"{program_id}_evaluable_genes_v1.tsv"
    )
    require(gene_file)

    eg = pd.read_csv(gene_file, sep="\t", dtype=str).fillna("")
    genes = [canon_symbol(x) for x in eg["Hugo_Symbol"]]

    tcols = np.array(
        [target["gene_index"][g] for g in genes],
        dtype=np.int64,
    )
    target_x = target["x"][:, tcols]
    target_z = target["z"][:, tcols]

    wdf = weights_df.loc[
        weights_df["program_id"] == program_id
    ].copy()
    wdf["Hugo_Symbol"] = wdf["Hugo_Symbol"].map(canon_symbol)
    wdf = wdf.set_index("Hugo_Symbol", drop=False)

    source_loading = np.array(
        [float(wdf.loc[g, "source_pc1_loading"]) for g in genes],
        dtype=np.float64,
    )
    source_sign = np.array(
        [float(wdf.loc[g, "source_pc1_loading_sign"]) for g in genes],
        dtype=np.float64,
    )

    source_loading_rank, source_loading_norm = centered_rank(source_loading)
    source_edge_r, source_edge_norm, tri_i, tri_j = source_edge_rank(
        source,
        genes,
    )

    x_t = torch.as_tensor(target_x, dtype=torch.float64, device=device)
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

    # Uniform-weight replay is mandatory.
    uniform = torch.full(
        (target_x.shape[0],),
        1.0 / target_x.shape[0],
        dtype=torch.float64,
        device=device,
    )
    corr_u, z_u = weighted_corr(x_t, uniform)
    edges_u = corr_u[ii_t, jj_t].detach().cpu().numpy()

    replay_edge = spearman_against_fixed_rank(
        source_edge_r,
        source_edge_norm,
        edges_u,
    )
    load_u_t, replay_resid, replay_fallback = leading_pc1_weighted(
        corr_u,
        z_u,
        uniform,
        init_t,
        sign_t,
    )
    load_u = load_u_t.detach().cpu().numpy()
    replay_load = spearman_against_fixed_rank(
        source_loading_rank,
        source_loading_norm,
        load_u,
    )

    observed_edge = float(row["rho_edge"])
    observed_load = float(row["rho_load"])

    if abs(replay_edge - observed_edge) > REPLAY_TOL:
        raise RuntimeError(
            f"{target_name} {program_id}: uniform-weight edge replay mismatch "
            f"{replay_edge} vs {observed_edge}"
        )
    if abs(replay_load - observed_load) > REPLAY_TOL:
        raise RuntimeError(
            f"{target_name} {program_id}: uniform-weight loading replay mismatch "
            f"{replay_load} vs {observed_load}"
        )

    del corr_u, z_u, edges_u, load_u_t, uniform
    torch.cuda.empty_cache()

    rel = reliability_diagnostics(
        target_z=target_z,
        signs=source_sign,
        repeats=SPLIT_HALF_REPEATS,
        seed=combo_seed(target_name, program_id, 401),
        device=device,
    )

    edge_rho = np.full(BAYESIAN_DRAWS, np.nan, dtype=np.float64)
    load_rho = np.full(BAYESIAN_DRAWS, np.nan, dtype=np.float64)
    pc1_resid = np.full(BAYESIAN_DRAWS, np.nan, dtype=np.float64)
    pc1_fallback = np.zeros(BAYESIAN_DRAWS, dtype=np.int8)

    seed = combo_seed(target_name, program_id, 501)
    completed = 0

    if paths["checkpoint"].exists():
        ck = np.load(paths["checkpoint"])
        if (
            int(ck["seed"]) == seed
            and int(ck["draws"]) == BAYESIAN_DRAWS
        ):
            completed = int(ck["completed"])
            edge_rho[:] = ck["edge_rho"]
            load_rho[:] = ck["load_rho"]
            pc1_resid[:] = ck["pc1_resid"]
            pc1_fallback[:] = ck["pc1_fallback"]
            print(f"    resuming Bayesian bootstrap from {completed}/{BAYESIAN_DRAWS}")

    rng = np.random.default_rng(seed)
    for _ in range(completed):
        _ = rng.exponential(scale=1.0, size=target_x.shape[0])

    t0 = time.perf_counter()

    for b in range(completed, BAYESIAN_DRAWS):
        w_np = rng.exponential(
            scale=1.0,
            size=target_x.shape[0],
        ).astype(np.float64)
        w_np /= w_np.sum()

        wt = torch.as_tensor(
            w_np,
            dtype=torch.float64,
            device=device,
        )

        corr, z = weighted_corr(x_t, wt)

        edges = corr[ii_t, jj_t].detach().cpu().numpy()
        edge_rho[b] = spearman_against_fixed_rank(
            source_edge_r,
            source_edge_norm,
            edges,
        )

        vb_t, resid, fallback = leading_pc1_weighted(
            corr,
            z,
            wt,
            init_t,
            sign_t,
        )
        vb = vb_t.detach().cpu().numpy()
        load_rho[b] = spearman_against_fixed_rank(
            source_loading_rank,
            source_loading_norm,
            vb,
        )
        pc1_resid[b] = resid
        pc1_fallback[b] = int(fallback)

        del wt, corr, z, vb_t

        if (
            (b + 1) % CHECKPOINT_EVERY == 0
            or (b + 1) == BAYESIAN_DRAWS
        ):
            np.savez_compressed(
                paths["checkpoint"],
                seed=np.array(seed, dtype=np.int64),
                draws=np.array(BAYESIAN_DRAWS, dtype=np.int64),
                completed=np.array(b + 1, dtype=np.int64),
                edge_rho=edge_rho,
                load_rho=load_rho,
                pc1_resid=pc1_resid,
                pc1_fallback=pc1_fallback,
            )
            print(
                f"    draws {b+1:4d}/{BAYESIAN_DRAWS}; "
                f"edge median={np.nanmedian(edge_rho[:b+1]):+.4f}; "
                f"load median={np.nanmedian(load_rho[:b+1]):+.4f}; "
                f"elapsed={time.perf_counter()-t0:.1f}s"
            )

    alpha = (1.0 - INTERVAL) / 2.0

    result = {
        "script_version": SCRIPT_VERSION,
        "target": target_name,
        "program_id": program_id,
        "n_samples": int(target_x.shape[0]),
        "n_evaluable_genes": len(genes),
        "n_edges": len(tri_i),
        "bayesian_bootstrap_draws": BAYESIAN_DRAWS,
        "bayesian_bootstrap_seed": seed,
        "bayesian_bootstrap_weight_distribution": "Dirichlet(1,...,1) via normalized Exp(1)",
        "observed_rho_edge": observed_edge,
        "observed_rho_load": observed_load,
        "uniform_weight_replay_rho_edge": replay_edge,
        "uniform_weight_replay_rho_load": replay_load,
        "uniform_weight_pc1_eigen_residual": replay_resid,
        "uniform_weight_pc1_exact_eigh_fallback": int(replay_fallback),
        "edge_bayesian_bootstrap_median": float(np.median(edge_rho)),
        "edge_bayesian_interval_low": float(np.quantile(edge_rho, alpha)),
        "edge_bayesian_interval_high": float(np.quantile(edge_rho, 1.0-alpha)),
        "loading_bayesian_bootstrap_median": float(np.median(load_rho)),
        "loading_bayesian_interval_low": float(np.quantile(load_rho, alpha)),
        "loading_bayesian_interval_high": float(np.quantile(load_rho, 1.0-alpha)),
        "loading_exact_eigh_fallbacks": int(pc1_fallback.sum()),
        "loading_max_eigen_residual": float(np.max(pc1_resid)),
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
        "leave_one_out_min": float(np.min(rel["leave_one_out_correlations"])),
        "leave_one_out_median": float(np.median(rel["leave_one_out_correlations"])),
        "elapsed_bayesian_bootstrap_seconds": time.perf_counter() - t0,
        "primary_classification_changed": False,
    }

    np.savez_compressed(
        paths["arrays"],
        edge_bayesian_bootstrap_rho=edge_rho,
        loading_bayesian_bootstrap_rho=load_rho,
        loading_pc1_eigen_residual=pc1_resid,
        loading_exact_eigh_fallback=pc1_fallback,
        split_half_raw=rel["split_half_raw"],
        split_half_spearman_brown=rel["split_half_spearman_brown"],
        leave_one_out_correlations=rel["leave_one_out_correlations"],
    )

    paths["result"].write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )

    del x_t, init_t, sign_t, ii_t, jj_t
    torch.cuda.empty_cache()

    return result


def main() -> None:
    print("=" * 142)
    print("Paper 4 / TCBB - Bayesian-bootstrap pooled target uncertainty + reliability diagnostics")
    print("=" * 142)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Frozen scientific contract:")
    print("  Target outcomes/treatment loaded:                    NO")
    print(f"  Bayesian-bootstrap draws per assessable combination: {BAYESIAN_DRAWS:,}")
    print("  Weight distribution:                                 Dirichlet(1,...,1)")
    print("  Source object:                                       FIXED")
    print("  Module genes:                                        FIXED")
    print("  Edge statistic:                                      full evaluable-edge rho_edge")
    print("  Loading statistic:                                   weighted target PC1 -> rho_load")
    print(f"  Uncertainty interval:                                {int(INTERVAL*100)}% Bayesian-bootstrap percentile")
    print(f"  Split-half repeats:                                  {SPLIT_HALF_REPEATS:,}")
    print("  Primary classification:                              FROZEN, NOT CHANGED")
    print("=" * 142)

    for p in [
        BAYESIAN_AMENDMENT,
        FINAL_CLASSIFICATION,
        UNIVERSE,
        WEIGHTS,
        TCGA_EXPR,
        SCANB_EXPR,
        SCANB_PRIMARY,
        METABRIC_EXPR,
    ]:
        require(p)

    acon = json.loads(BAYESIAN_AMENDMENT.read_text(encoding="utf-8"))
    if acon.get("status") != "FROZEN_BEFORE_ANY_VALID_TARGET_BOOTSTRAP_PRESERVATION_RESULT":
        raise RuntimeError("05c2 Bayesian-bootstrap amendment has unexpected status.")
    if int(acon["bayesian_bootstrap"]["draws"]) != BAYESIAN_DRAWS:
        raise RuntimeError("Bayesian-bootstrap draw count differs from frozen amendment.")

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
        raise RuntimeError(f"Expected 10,000 frozen source genes, got {len(universe)}.")

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

    print("\n[1/4] Replaying frozen TCGA source ...")
    source = load_tcga_source(universe)
    print(f"  {source['z'].shape[0]:,} samples x {source['z'].shape[1]:,} genes")

    results = []

    for target_name in ["SCANB_GSE96058", "METABRIC"]:
        print(f"\n[2/4] Loading {target_name} ...")
        target = (
            load_scanb_target(universe_set, primary_titles)
            if target_name == "SCANB_GSE96058"
            else load_metabric_target(universe_set)
        )
        print(
            f"  {target['x'].shape[0]:,} samples x "
            f"{target['x'].shape[1]:,} evaluable genes"
        )

        print("\n" + "=" * 142)
        print(f"[3/4] BAYESIAN BOOTSTRAP + RELIABILITY — {target_name}")
        print("=" * 142)

        sub = classification.loc[
            classification["target"] == target_name
        ].copy()

        for i, row in enumerate(sub.to_dict(orient="records"), start=1):
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
                        "reason": "frozen primary assessability guard failed",
                    }
                )
                continue

            print(
                f"\n  [{i:02d}/12] {program_id}: "
                f"CLASS={row['primary_classification']}; "
                f"observed edge={float(row['rho_edge']):+.4f}; "
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
                f"    edge 95% Bayesian interval="
                f"[{result['edge_bayesian_interval_low']:+.4f},"
                f"{result['edge_bayesian_interval_high']:+.4f}]; "
                f"loading="
                f"[{result['loading_bayesian_interval_low']:+.4f},"
                f"{result['loading_bayesian_interval_high']:+.4f}]"
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

    print("\n[4/4] Writing Bayesian-bootstrap/reliability summary ...")

    df = pd.DataFrame(results)
    summary_tsv = OUT_DIR / "primary_pooled_bayesian_bootstrap_reliability_summary_v3.tsv"
    df.to_csv(summary_tsv, sep="\t", index=False)

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "PRIMARY_POOLED_BAYESIAN_BOOTSTRAP_RELIABILITY_COMPLETE",
        "bayesian_bootstrap_draws": BAYESIAN_DRAWS,
        "interval": INTERVAL,
        "split_half_repeats": SPLIT_HALF_REPEATS,
        "ordinary_bootstrap_failure_retained_as_diagnostic": True,
        "primary_classification_changed": False,
        "summary_file": str(summary_tsv),
    }
    master_json = OUT_DIR / "primary_pooled_bayesian_bootstrap_reliability_v3.json"
    master_json.write_text(json.dumps(master, indent=2), encoding="utf-8")

    print("\n" + "=" * 142)
    print("05c v3 BAYESIAN-BOOTSTRAP / RELIABILITY: COMPLETE")
    print("=" * 142)
    print("Primary classifications were NOT changed.")
    print("Ordinary fixed-gene bootstrap failure remains documented as a diagnostic.")
    print("No target outcomes or treatment variables were loaded.")
    print()
    print(f"Summary TSV: {summary_tsv}")
    print(f"Master JSON: {master_json}")
    print("=" * 142)


if __name__ == "__main__":
    main()
