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
    raise RuntimeError("PyTorch is required for 05c. Run the CUDA audit first.") from exc

try:
    from scipy.stats import rankdata
except Exception as exc:
    raise RuntimeError("scipy is required for exact Spearman statistics.") from exc


SCRIPT_VERSION = "05c-run-primary-bootstrap-reliability-gpu-v2-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

PRESERVATION_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_preservation_operating_contract_v1"
    / "preservation_operating_characteristics_contract_v1.json"
)
BOOTSTRAP_DEGENERACY_AMENDMENT = (
    DATA_ROOT
    / "paper4_tcbb_bootstrap_degeneracy_amendment_v1"
    / "bootstrap_degeneracy_amendment_v1.json"
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

OUT_DIR = DATA_ROOT / "paper4_tcbb_primary_bootstrap_reliability_v2"

BOOTSTRAP_REPLICATES = 1000
MAX_BOOTSTRAP_ATTEMPTS = 10000
BOOTSTRAP_CI = 0.95
SPLIT_HALF_REPEATS = 2000
BASE_SEED = 20260916

CHECKPOINT_EVERY = 25
SPLIT_BATCH = 200

PC1_RESIDUAL_TOL = 1e-10
MAX_POWER_ITER = 400
OBSERVED_REPLAY_TOL = 5e-10

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


def centered_rank(x: np.ndarray) -> tuple[np.ndarray, float]:
    r = rankdata(np.asarray(x, dtype=np.float64), method="average")
    r -= r.mean()
    ss = float(np.dot(r, r))
    if ss <= 0:
        raise RuntimeError("Degenerate rank vector.")
    return r, math.sqrt(ss)


def load_tcga_source(universe: pd.DataFrame) -> dict:
    print("\n[1/4] Replaying frozen TCGA source matrix ...")
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
    print("\n[2/4] Loading frozen SCAN-B primary target matrix ...")
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
    x = x[:, eligible]

    if len(genes) != 9225:
        raise RuntimeError(f"Expected 9,225 SCAN-B evaluable genes, got {len(genes)}.")

    print(f"  {x.shape[0]:,} samples x {x.shape[1]:,} evaluable genes; {time.perf_counter()-t0:.1f}s")
    return {
        "name": "SCANB_GSE96058",
        "genes": genes,
        "gene_index": {g: i for i, g in enumerate(genes)},
        "x": x,
        "z": standardize_samples_x_genes_np(x),
    }


def load_metabric_target(universe_set: set[str]) -> dict:
    print("\n[2/4] Loading frozen METABRIC target matrix ...")
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
    x = x[:, eligible]

    if len(genes) != 8485:
        raise RuntimeError(f"Expected 8,485 METABRIC evaluable genes, got {len(genes)}.")

    print(f"  {x.shape[0]:,} samples x {x.shape[1]:,} evaluable genes; {time.perf_counter()-t0:.1f}s")
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
    cols = np.array([source["gene_index"][g] for g in genes], dtype=np.int64)
    z = source["z"][:, cols]
    corr = (z.T @ z) / (z.shape[0] - 1)
    tri_i, tri_j = np.triu_indices(len(genes), k=1)
    edges = corr[tri_i, tri_j]
    rank = rankdata(edges, method="average")
    rank -= rank.mean()
    norm = math.sqrt(float(np.dot(rank, rank)))
    if norm <= 0:
        raise RuntimeError("Degenerate frozen source edge-rank vector.")
    return rank, norm, tri_i, tri_j


def leading_pc1_from_corr(
    corr_t: torch.Tensor,
    z_t: torch.Tensor,
    init_loading: torch.Tensor,
    sign_t: torch.Tensor,
) -> tuple[torch.Tensor, float, bool]:
    eps = torch.finfo(torch.float64).eps
    v = init_loading / torch.linalg.vector_norm(init_loading)

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

    scores_c = scores - scores.mean()
    sign_c = sign_score - sign_score.mean()
    den = torch.linalg.vector_norm(scores_c) * torch.linalg.vector_norm(sign_c)
    if den <= 0:
        raise RuntimeError("Degenerate PC1 orientation score.")
    orient = torch.dot(scores_c, sign_c) / den
    if orient < 0:
        v = -v

    return v, residual_value, fallback


def reliability_diagnostics(
    target_z: np.ndarray,
    signs: np.ndarray,
    repeats: int,
    seed: int,
    device: torch.device,
) -> dict:
    n, m = target_z.shape
    if m < 2:
        raise RuntimeError("Need at least two genes for split-half reliability.")

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

    # Leave-one-gene-out frozen-sign score correlation with the full score.
    full = torch.mean(signed, dim=1)
    loo = (
        m * full[:, None] - signed
    ) / (m - 1)

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


def checkpoint_paths(target: str, program_id: str) -> dict[str, Path]:
    d = OUT_DIR / "per_program" / target
    d.mkdir(parents=True, exist_ok=True)
    return {
        "checkpoint": d / f"{program_id}_bootstrap_checkpoint_v2.npz",
        "final_npz": d / f"{program_id}_bootstrap_reliability_arrays_v2.npz",
        "result": d / f"{program_id}_bootstrap_reliability_result_v2.json",
        "degenerate": d / f"{program_id}_bootstrap_degenerate_gene_counts_v2.tsv",
    }


def bootstrap_program(
    target: dict,
    source: dict,
    row: pd.Series,
    weights: pd.DataFrame,
    device: torch.device,
) -> dict:
    target_name = str(row["target"])
    program_id = str(row["program_id"])
    paths = checkpoint_paths(target_name, program_id)

    if paths["result"].exists() and paths["final_npz"].exists():
        try:
            old = json.loads(paths["result"].read_text(encoding="utf-8"))
        except Exception:
            old = {}
        if (
            old.get("script_version") == SCRIPT_VERSION
            and int(old.get("bootstrap_replicates", -1)) == BOOTSTRAP_REPLICATES
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
    m = len(genes)

    tcols = np.array([target["gene_index"][g] for g in genes], dtype=np.int64)
    target_x = target["x"][:, tcols]
    target_z_observed = target["z"][:, tcols]

    wdf = weights.loc[weights["program_id"] == program_id].copy()
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

    source_loading_rank = rankdata(source_loading, method="average")
    source_loading_rank -= source_loading_rank.mean()
    source_loading_rank_norm = math.sqrt(
        float(np.dot(source_loading_rank, source_loading_rank))
    )

    source_edge_rank, source_edge_rank_norm, tri_i, tri_j = (
        source_edge_rank_for_program(source, genes)
    )

    observed_edge = float(row["rho_edge"])
    observed_load = float(row["rho_load"])

    # Replay observed target edge/load statistics from the same target matrix before bootstrap.
    z_obs_t = torch.as_tensor(
        target_z_observed,
        dtype=torch.float64,
        device=device,
    )
    corr_obs_t = (z_obs_t.T @ z_obs_t) / (z_obs_t.shape[0] - 1)
    edges_obs = corr_obs_t[
        torch.as_tensor(tri_i, dtype=torch.long, device=device),
        torch.as_tensor(tri_j, dtype=torch.long, device=device),
    ].detach().cpu().numpy()

    replay_edge = spearman_against_fixed_rank(
        source_edge_rank,
        source_edge_rank_norm,
        edges_obs,
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
    loading_obs_t, obs_resid, obs_fallback = leading_pc1_from_corr(
        corr_obs_t,
        z_obs_t,
        init_t,
        sign_t,
    )
    loading_obs = loading_obs_t.detach().cpu().numpy()
    replay_load = spearman_against_fixed_rank(
        source_loading_rank,
        source_loading_rank_norm,
        loading_obs,
    )

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

    del corr_obs_t, edges_obs, loading_obs_t, z_obs_t
    torch.cuda.empty_cache()

    # Reliability is based on the original pooled standardized target matrix.
    rel_seed = combo_seed(target_name, program_id, 401)
    rel = reliability_diagnostics(
        target_z=target_z_observed,
        signs=source_sign,
        repeats=SPLIT_HALF_REPEATS,
        seed=rel_seed,
        device=device,
    )

    # Bootstrap checkpoint state. Arrays index VALID bootstrap draws, not attempts.
    edge_rho = np.full(BOOTSTRAP_REPLICATES, np.nan, dtype=np.float64)
    load_rho = np.full(BOOTSTRAP_REPLICATES, np.nan, dtype=np.float64)
    pc1_resid = np.full(BOOTSTRAP_REPLICATES, np.nan, dtype=np.float64)
    pc1_fallback = np.zeros(BOOTSTRAP_REPLICATES, dtype=np.int8)

    completed = 0
    attempts = 0
    invalid_attempts = 0
    degenerate_gene_counts = np.zeros(m, dtype=np.int64)

    bootstrap_seed = combo_seed(target_name, program_id, 301)

    if paths["checkpoint"].exists():
        ck = np.load(paths["checkpoint"])
        if (
            int(ck["bootstrap_seed"]) == bootstrap_seed
            and int(ck["bootstrap_replicates"]) == BOOTSTRAP_REPLICATES
            and int(ck["maximum_attempts"]) == MAX_BOOTSTRAP_ATTEMPTS
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
                f"    resuming bootstrap from valid={completed}/{BOOTSTRAP_REPLICATES}, "
                f"attempts={attempts}, invalid={invalid_attempts}"
            )

    rng = np.random.default_rng(bootstrap_seed)

    # Replay RNG state deterministically by ATTEMPT count.
    for _ in range(attempts):
        _ = rng.integers(
            0,
            target_x.shape[0],
            size=target_x.shape[0],
            endpoint=False,
        )

    x_full_t = torch.as_tensor(
        target_x,
        dtype=torch.float64,
        device=device,
    )
    tri_i_t = torch.as_tensor(tri_i, dtype=torch.long, device=device)
    tri_j_t = torch.as_tensor(tri_j, dtype=torch.long, device=device)

    t0 = time.perf_counter()
    last_checkpoint_completed = completed

    while (
        completed < BOOTSTRAP_REPLICATES
        and attempts < MAX_BOOTSTRAP_ATTEMPTS
    ):
        idx_np = rng.integers(
            0,
            target_x.shape[0],
            size=target_x.shape[0],
            endpoint=False,
        )
        attempts += 1

        idx_t = torch.as_tensor(
            idx_np,
            dtype=torch.long,
            device=device,
        )
        xb = x_full_t.index_select(0, idx_t)

        mu = xb.mean(dim=0, keepdim=True)
        sd = xb.std(dim=0, correction=1, keepdim=True)

        bad_t = (~torch.isfinite(sd)) | (sd <= 0)
        if bool(torch.any(bad_t)):
            bad = bad_t.squeeze(0).detach().cpu().numpy().astype(bool)
            degenerate_gene_counts[bad] += 1
            invalid_attempts += 1

            del idx_t, xb, mu, sd, bad_t
            continue

        zb = (xb - mu) / sd
        corr = (zb.T @ zb) / (zb.shape[0] - 1)

        edges_t = corr[tri_i_t, tri_j_t]
        edges = edges_t.detach().cpu().numpy()
        erank = rankdata(edges, method="average")
        erank -= erank.mean()
        enorm = math.sqrt(float(np.dot(erank, erank)))
        edge_rho[completed] = float(
            np.dot(source_edge_rank, erank)
            / (source_edge_rank_norm * enorm)
        )

        vb_t, resid, fallback = leading_pc1_from_corr(
            corr,
            zb,
            init_t,
            sign_t,
        )
        vb = vb_t.detach().cpu().numpy()
        vrank = rankdata(vb, method="average")
        vrank -= vrank.mean()
        vnorm = math.sqrt(float(np.dot(vrank, vrank)))
        load_rho[completed] = float(
            np.dot(source_loading_rank, vrank)
            / (source_loading_rank_norm * vnorm)
        )
        pc1_resid[completed] = resid
        pc1_fallback[completed] = int(fallback)

        completed += 1

        del (
            idx_t,
            xb,
            mu,
            sd,
            bad_t,
            zb,
            corr,
            edges_t,
            vb_t,
        )

        if (
            completed - last_checkpoint_completed >= CHECKPOINT_EVERY
            or completed == BOOTSTRAP_REPLICATES
        ):
            np.savez_compressed(
                paths["checkpoint"],
                bootstrap_seed=np.array(bootstrap_seed, dtype=np.int64),
                bootstrap_replicates=np.array(
                    BOOTSTRAP_REPLICATES,
                    dtype=np.int64,
                ),
                maximum_attempts=np.array(
                    MAX_BOOTSTRAP_ATTEMPTS,
                    dtype=np.int64,
                ),
                completed=np.array(completed, dtype=np.int64),
                attempts=np.array(attempts, dtype=np.int64),
                invalid_attempts=np.array(
                    invalid_attempts,
                    dtype=np.int64,
                ),
                edge_rho=edge_rho,
                load_rho=load_rho,
                pc1_resid=pc1_resid,
                pc1_fallback=pc1_fallback,
                degenerate_gene_counts=degenerate_gene_counts,
            )
            last_checkpoint_completed = completed

            elapsed = time.perf_counter() - t0
            print(
                f"    valid bootstraps {completed:4d}/{BOOTSTRAP_REPLICATES}; "
                f"attempts={attempts:4d}; invalid={invalid_attempts:4d}; "
                f"edge median={np.nanmedian(edge_rho[:completed]):+.4f}; "
                f"load median={np.nanmedian(load_rho[:completed]):+.4f}; "
                f"elapsed={elapsed:.1f}s"
            )

    if completed != BOOTSTRAP_REPLICATES:
        # Persist the failure diagnostics before stopping.
        np.savez_compressed(
            paths["checkpoint"],
            bootstrap_seed=np.array(bootstrap_seed, dtype=np.int64),
            bootstrap_replicates=np.array(
                BOOTSTRAP_REPLICATES,
                dtype=np.int64,
            ),
            maximum_attempts=np.array(
                MAX_BOOTSTRAP_ATTEMPTS,
                dtype=np.int64,
            ),
            completed=np.array(completed, dtype=np.int64),
            attempts=np.array(attempts, dtype=np.int64),
            invalid_attempts=np.array(
                invalid_attempts,
                dtype=np.int64,
            ),
            edge_rho=edge_rho,
            load_rho=load_rho,
            pc1_resid=pc1_resid,
            pc1_fallback=pc1_fallback,
            degenerate_gene_counts=degenerate_gene_counts,
        )

        pd.DataFrame(
            {
                "Hugo_Symbol": genes,
                "degenerate_attempt_count": degenerate_gene_counts,
            }
        ).sort_values(
            "degenerate_attempt_count",
            ascending=False,
        ).to_csv(
            paths["degenerate"],
            sep="\t",
            index=False,
        )

        raise RuntimeError(
            f"{target_name} {program_id}: only {completed} valid bootstrap "
            f"draws were obtained in {attempts} attempts. Frozen cap is "
            f"{MAX_BOOTSTRAP_ATTEMPTS}; uncertainty is not estimable under "
            f"the frozen nonparametric scheme."
        )

    pd.DataFrame(
        {
            "Hugo_Symbol": genes,
            "degenerate_attempt_count": degenerate_gene_counts,
        }
    ).sort_values(
        "degenerate_attempt_count",
        ascending=False,
    ).to_csv(
        paths["degenerate"],
        sep="\t",
        index=False,
    )

    del (
        x_full_t,
        tri_i_t,
        tri_j_t,
        init_t,
        sign_t,
    )
    torch.cuda.empty_cache()

    alpha = (1.0 - BOOTSTRAP_CI) / 2.0

    result = {
        "script_version": SCRIPT_VERSION,
        "target": target_name,
        "program_id": program_id,
        "n_samples": int(target_x.shape[0]),
        "n_evaluable_genes": m,
        "n_edges": int(len(tri_i)),
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "bootstrap_seed": bootstrap_seed,
        "bootstrap_attempts": attempts,
        "bootstrap_invalid_attempts": invalid_attempts,
        "bootstrap_valid_attempt_fraction": float(
            BOOTSTRAP_REPLICATES / attempts
        ),
        "maximum_bootstrap_attempts": MAX_BOOTSTRAP_ATTEMPTS,
        "genes_ever_degenerate_in_bootstrap": int(
            np.sum(degenerate_gene_counts > 0)
        ),
        "maximum_single_gene_degenerate_attempt_count": int(
            np.max(degenerate_gene_counts)
        ),
        "observed_rho_edge": observed_edge,
        "observed_rho_load": observed_load,
        "observed_replay_rho_edge": replay_edge,
        "observed_replay_rho_load": replay_load,
        "observed_pc1_eigen_residual": obs_resid,
        "observed_pc1_exact_eigh_fallback": int(obs_fallback),
        "edge_bootstrap_median": float(np.median(edge_rho)),
        "edge_bootstrap_ci_low": float(np.quantile(edge_rho, alpha)),
        "edge_bootstrap_ci_high": float(np.quantile(edge_rho, 1.0 - alpha)),
        "loading_bootstrap_median": float(np.median(load_rho)),
        "loading_bootstrap_ci_low": float(np.quantile(load_rho, alpha)),
        "loading_bootstrap_ci_high": float(np.quantile(load_rho, 1.0 - alpha)),
        "loading_bootstrap_exact_eigh_fallbacks": int(pc1_fallback.sum()),
        "loading_bootstrap_max_eigen_residual": float(np.max(pc1_resid)),
        "split_half_repeats": SPLIT_HALF_REPEATS,
        "split_half_seed": rel_seed,
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
        paths["final_npz"],
        edge_bootstrap_rho=edge_rho,
        loading_bootstrap_rho=load_rho,
        loading_pc1_eigen_residual=pc1_resid,
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

    return result


def main() -> None:
    print("=" * 142)
    print("Paper 4 / TCBB - primary pooled target-sampling bootstrap uncertainty + reliability diagnostics")
    print("=" * 142)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Frozen scientific contract:")
    print("  Target outcomes/treatment loaded:                    NO")
    print(f"  VALID target-sample bootstraps per combination:      {BOOTSTRAP_REPLICATES:,}")
    print(f"  Maximum bootstrap attempts per combination:          {MAX_BOOTSTRAP_ATTEMPTS:,}")
    print("  Degenerate resample handling:                        discard undefined attempt + redraw")
    print("  Per-replicate gene dropping/jitter/epsilon:           PROHIBITED")
    print("  Source object during target bootstrap:               FIXED")
    print("  Edge bootstrap statistic:                            full evaluable-edge rho_edge")
    print("  Loading bootstrap:                                   recompute target PC1 + rho_load")
    print(f"  Percentile CI:                                       {int(BOOTSTRAP_CI*100)}%")
    print(f"  Split-half repeats:                                  {SPLIT_HALF_REPEATS:,}")
    print("  Split-half classifier threshold:                     NONE")
    print("  Leave-one-out:                                       descriptive")
    print("  Primary Strong/Partial/No-clear classification:      FROZEN, NOT CHANGED")
    print("=" * 142)

    for p in [
        PRESERVATION_CONTRACT,
        BOOTSTRAP_DEGENERACY_AMENDMENT,
        FINAL_CLASSIFICATION,
        UNIVERSE,
        WEIGHTS,
        TCGA_EXPR,
        SCANB_EXPR,
        SCANB_PRIMARY,
        METABRIC_EXPR,
    ]:
        require(p)

    contract = json.loads(
        PRESERVATION_CONTRACT.read_text(encoding="utf-8")
    )
    degeneracy_amendment = json.loads(
        BOOTSTRAP_DEGENERACY_AMENDMENT.read_text(encoding="utf-8")
    )
    if contract.get("status") != "FROZEN_BEFORE_FIRST_TARGET_PRESERVATION_STATISTIC":
        raise RuntimeError("04a preservation contract has unexpected status.")

    if int(contract["target_sampling_uncertainty"]["bootstrap_replicates"]) != BOOTSTRAP_REPLICATES:
        raise RuntimeError("Bootstrap replicate count differs from frozen 04a.")
    if int(contract["reliability"]["split_half_repeats"]) != SPLIT_HALF_REPEATS:
        raise RuntimeError("Split-half repeat count differs from frozen 04a.")

    if degeneracy_amendment.get("status") != "FROZEN_BEFORE_ANY_BOOTSTRAP_PRESERVATION_RESULT":
        raise RuntimeError("05c0 degeneracy amendment has unexpected status.")
    if int(
        degeneracy_amendment["degenerate_draw_rule"]["required_valid_draws"]
    ) != BOOTSTRAP_REPLICATES:
        raise RuntimeError("05c0 valid-draw requirement differs from runner.")
    if int(
        degeneracy_amendment["degenerate_draw_rule"]["maximum_total_attempts"]
    ) != MAX_BOOTSTRAP_ATTEMPTS:
        raise RuntimeError("05c0 attempt cap differs from runner.")

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
            f"Expected 10,000 frozen source genes, got {len(universe)}."
        )

    universe_genes = [
        canon_symbol(x)
        for x in universe["Hugo_Symbol"]
    ]
    universe_set = set(universe_genes)

    weights = pd.read_csv(
        WEIGHTS,
        sep="\t",
        dtype=str,
    ).fillna("")
    weights["Hugo_Symbol"] = weights["Hugo_Symbol"].map(
        canon_symbol
    )

    primary = pd.read_csv(
        SCANB_PRIMARY,
        sep="\t",
        dtype=str,
    ).fillna("")
    primary_titles = list(primary["primary_title"])

    source = load_tcga_source(universe)

    results = []

    for target_name in [
        "SCANB_GSE96058",
        "METABRIC",
    ]:
        target = (
            load_scanb_target(
                universe_set,
                primary_titles,
            )
            if target_name == "SCANB_GSE96058"
            else load_metabric_target(universe_set)
        )

        print("\n" + "=" * 142)
        print(
            f"[3/4] BOOTSTRAP + RELIABILITY — {target_name}"
        )
        print("=" * 142)

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
                    "NOT ASSESSABLE — bootstrap/reliability not estimated"
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

            result = bootstrap_program(
                target=target,
                source=source,
                row=pd.Series(row),
                weights=weights,
                device=device,
            )
            results.append(result)

            print(
                f"    edge {int(BOOTSTRAP_CI*100)}% CI="
                f"[{result['edge_bootstrap_ci_low']:+.4f},"
                f"{result['edge_bootstrap_ci_high']:+.4f}]; "
                f"loading CI="
                f"[{result['loading_bootstrap_ci_low']:+.4f},"
                f"{result['loading_bootstrap_ci_high']:+.4f}]"
            )
            print(
                f"    split-half raw median="
                f"{result['split_half_raw_median']:.4f}; "
                f"Spearman-Brown median="
                f"{result['split_half_spearman_brown_median']:.4f}; "
                f"LOO min={result['leave_one_out_min']:.6f}"
            )

        del target
        torch.cuda.empty_cache()
        gc.collect()

    print("\n[4/4] Writing pooled bootstrap/reliability summary ...")

    out_df = pd.DataFrame(results)
    summary_tsv = OUT_DIR / "primary_pooled_bootstrap_reliability_summary_v2.tsv"
    out_df.to_csv(
        summary_tsv,
        sep="\t",
        index=False,
    )

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "PRIMARY_POOLED_BOOTSTRAP_RELIABILITY_COMPLETE_V2",
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "maximum_bootstrap_attempts": MAX_BOOTSTRAP_ATTEMPTS,
        "bootstrap_degeneracy_amendment": str(BOOTSTRAP_DEGENERACY_AMENDMENT),
        "bootstrap_ci": BOOTSTRAP_CI,
        "split_half_repeats": SPLIT_HALF_REPEATS,
        "primary_classification_changed": False,
        "summary_file": str(summary_tsv),
    }

    master_json = OUT_DIR / "primary_pooled_bootstrap_reliability_v2.json"
    master_json.write_text(
        json.dumps(master, indent=2),
        encoding="utf-8",
    )

    print("\n" + "=" * 142)
    print("05c v2 PRIMARY POOLED BOOTSTRAP / RELIABILITY: COMPLETE")
    print("=" * 142)
    print("Primary classifications were NOT changed.")
    print("No target outcomes or treatment variables were loaded.")
    print()
    print(f"Summary TSV: {summary_tsv}")
    print(f"Master JSON: {master_json}")
    print("=" * 142)


if __name__ == "__main__":
    main()
