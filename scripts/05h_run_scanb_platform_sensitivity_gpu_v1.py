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
    raise RuntimeError("PyTorch is required for 05h.") from exc

try:
    from scipy.stats import rankdata
except Exception as exc:
    raise RuntimeError("scipy is required for Spearman statistics.") from exc


SCRIPT_VERSION = "05h-run-scanb-platform-sensitivity-gpu-v1-no-cli"

DATA_ROOT = Path(r"D:\paper4_tcbb_data")

EXEC_CONTRACT = (
    DATA_ROOT
    / "paper4_tcbb_scanb_platform_sensitivity_execution_contract_v1"
    / "scanb_platform_sensitivity_execution_contract_v1.json"
)
PLATFORM_MANIFEST = (
    DATA_ROOT
    / "paper4_tcbb_pam50_alignment_audit_v1"
    / "scanb_frozen_primary_pam50_alignment_v1.tsv"
)
PLATFORM_COMPOSITION = (
    DATA_ROOT
    / "paper4_tcbb_final_pretarget_audits_v1"
    / "scanb_platform_pam50_composition_v1.tsv"
)
CLASSIFICATION = (
    DATA_ROOT
    / "paper4_tcbb_final_mapping_specificity_null_v2"
    / "primary_pooled_final_classification_v2.tsv"
)
SPECIFICITY_ROOT = DATA_ROOT / "paper4_tcbb_final_mapping_specificity_null_v2"
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

OUT_DIR = DATA_ROOT / "paper4_tcbb_scanb_platform_sensitivity_v1"

PLATFORM_COUNTS = {
    "GPL11154": 2969,
    "GPL18573": 304,
}
PLATFORM_CODE = {
    "GPL11154": 1,
    "GPL18573": 2,
}
BASE_SEED = 20260918
BOOTSTRAP_REPLICATES = 1000
MAX_ATTEMPTS = 10000
CHECKPOINT_EVERY = 50

EDGE_BATCH = 1024
CPU_VALIDATION_EDGES = 25
CPU_GPU_TOL = 3e-10
PC1_RESIDUAL_TOL = 1e-10
MAX_POWER_ITER = 400


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


def centered_rank_cpu(x: np.ndarray) -> tuple[np.ndarray, float]:
    r = rankdata(np.asarray(x, dtype=np.float64), method="average")
    r -= r.mean()
    ss = float(np.dot(r, r))
    if ss <= 0:
        raise RuntimeError("Degenerate rank vector.")
    return r, math.sqrt(ss)


def spearman_cpu(x: np.ndarray, y: np.ndarray) -> float:
    xr, xn = centered_rank_cpu(x)
    yr, yn = centered_rank_cpu(y)
    return float(np.dot(xr, yr) / (xn * yn))


def gpu_centered_average_ranks(values: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    sorted_vals, order = torch.sort(values)
    n = sorted_vals.numel()

    new_group = torch.empty(n, dtype=torch.bool, device=values.device)
    new_group[0] = True
    if n > 1:
        new_group[1:] = sorted_vals[1:] != sorted_vals[:-1]

    group_id = torch.cumsum(new_group.to(torch.int64), dim=0) - 1
    counts = torch.bincount(group_id)
    starts0 = torch.cumsum(counts, dim=0) - counts
    avg_rank = starts0.to(torch.float64) + (counts.to(torch.float64) - 1.0) / 2.0 + 1.0
    ranks_sorted = avg_rank[group_id]

    ranks = torch.empty(n, dtype=torch.float64, device=values.device)
    ranks[order] = ranks_sorted
    centered = ranks - (n + 1.0) / 2.0
    norm = torch.linalg.vector_norm(centered)

    del sorted_vals, order, new_group, group_id, counts, starts0, avg_rank, ranks_sorted, ranks
    return centered, norm


def gpu_spearman_against_fixed(
    fixed_rank_t: torch.Tensor,
    fixed_norm_t: torch.Tensor,
    values_t: torch.Tensor,
) -> float:
    yr, yn = gpu_centered_average_ranks(values_t)
    rho = torch.dot(fixed_rank_t, yr) / (fixed_norm_t * yn)
    out = float(rho.detach().cpu())
    del yr, yn, rho
    return out


def standardize_np(x: np.ndarray) -> np.ndarray:
    if not np.isfinite(x).all():
        raise RuntimeError("Non-finite values reached standardization.")
    mu = x.mean(axis=0, keepdims=True)
    sd = x.std(axis=0, ddof=1, keepdims=True)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise RuntimeError("Zero/non-finite SD reached standardization.")
    return (x - mu) / sd


def exact_variable_np(x: np.ndarray) -> np.ndarray:
    finite = np.isfinite(x).all(axis=0)
    out = np.zeros(x.shape[1], dtype=bool)
    if finite.any():
        xmin = np.min(x[:, finite], axis=0)
        xmax = np.max(x[:, finite], axis=0)
        out[finite] = xmax > xmin
    return out


def full_corr_from_x(x: np.ndarray, device: torch.device) -> tuple[np.ndarray, float]:
    z = standardize_np(x)
    zt = torch.as_tensor(z, dtype=torch.float64, device=device)
    corr_t = (zt.T @ zt) / (z.shape[0] - 1)
    corr = corr_t.detach().cpu().numpy()

    k = min(30, x.shape[1])
    cpu = (z[:, :k].T @ z[:, :k]) / (z.shape[0] - 1)
    diff = float(np.max(np.abs(cpu - corr[:k, :k])))
    if diff > CPU_GPU_TOL:
        raise RuntimeError(
            f"CPU/GPU full-correlation validation failed: max|Δ|={diff:.3e}"
        )

    np.fill_diagonal(corr, 1.0)
    del zt, corr_t
    torch.cuda.empty_cache()
    return corr, diff


def pc1_loading(
    corr: np.ndarray,
    x: np.ndarray,
    source_loading: np.ndarray,
    source_sign: np.ndarray,
    random_start: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, dict]:
    z = standardize_np(x)
    ct = torch.as_tensor(corr, dtype=torch.float64, device=device)
    eps = torch.finfo(torch.float64).eps

    starts = [
        np.asarray(source_loading, dtype=np.float64),
        np.asarray(source_sign, dtype=np.float64),
        np.asarray(random_start, dtype=np.float64),
    ]
    candidates = []

    for start_idx, start in enumerate(starts):
        if not np.isfinite(start).all() or np.linalg.norm(start) == 0:
            continue

        v = torch.as_tensor(start, dtype=torch.float64, device=device)
        v = v / torch.linalg.vector_norm(v)
        residual_value = float("inf")

        for it in range(1, MAX_POWER_ITER + 1):
            w = ct @ v
            wn = torch.linalg.vector_norm(w)
            if not torch.isfinite(wn) or wn <= 0:
                break

            nv = w / wn
            if torch.dot(nv, v) < 0:
                nv = -nv
            v = nv

            if it % 5 == 0 or it == MAX_POWER_ITER:
                av = ct @ v
                lam = torch.dot(v, av)
                residual = torch.linalg.vector_norm(
                    av - lam * v
                ) / (torch.abs(lam) + eps)
                residual_value = float(residual.detach().cpu())
                if residual_value <= PC1_RESIDUAL_TOL:
                    break

        av = ct @ v
        lam = float(torch.dot(v, av).detach().cpu())
        candidates.append(
            (
                lam,
                residual_value,
                it,
                start_idx,
                v.detach().cpu().numpy(),
            )
        )

    if not candidates:
        raise RuntimeError("No valid target-PC1 power-iteration candidate.")

    candidates.sort(key=lambda q: q[0], reverse=True)
    lam, residual, iters, start_idx, loading = candidates[0]
    exact_fallback = False

    if residual > PC1_RESIDUAL_TOL:
        evals, evecs = torch.linalg.eigh(ct)
        v = evecs[:, -1]
        lam_t = evals[-1]
        av = ct @ v
        residual_t = torch.linalg.vector_norm(
            av - lam_t * v
        ) / (torch.abs(lam_t) + eps)
        loading = v.detach().cpu().numpy()
        lam = float(lam_t.detach().cpu())
        residual = float(residual_t.detach().cpu())
        exact_fallback = True

    scores = z @ loading
    sign_score = np.mean(z * source_sign[None, :], axis=1)
    orient = float(np.corrcoef(scores, sign_score)[0, 1])
    if not np.isfinite(orient):
        raise RuntimeError("Target PC1 orientation correlation is non-finite.")
    if orient < 0:
        loading = -loading
        orient = -orient

    meta = {
        "pc1_eigen_residual": residual,
        "pc1_exact_eigh_fallback": int(exact_fallback),
        "pc1_power_start_index": int(start_idx),
        "pc1_power_iterations": int(iters),
        "pc1_orientation_correlation_with_frozen_sign_score": orient,
        "pc1_variance_explained": float(lam / np.trace(corr)),
    }

    del ct
    torch.cuda.empty_cache()
    return loading, meta


def load_source_raw(universe: pd.DataFrame) -> dict:
    raw = pd.read_csv(TCGA_EXPR, sep="\t", low_memory=False)

    row_idx = universe["source_row_index_0based"].astype(int).to_numpy()
    selected = raw.iloc[row_idx, :]

    expected = [canon_symbol(x) for x in universe["Hugo_Symbol"]]
    observed = [canon_symbol(x) for x in selected["Hugo_Symbol"]]
    if expected != observed:
        raise RuntimeError("TCGA frozen source-universe replay mismatch.")

    values = selected.iloc[:, 2:].to_numpy(dtype=np.float64)
    if np.any(values < 0) or not np.isfinite(values).all():
        raise RuntimeError("TCGA RSEM contains invalid values.")

    x = np.log2(values + 1.0).T
    if x.shape[0] != 1082 or x.shape[1] != 10000:
        raise RuntimeError(f"Unexpected TCGA source shape: {x.shape}")

    return {
        "samples": list(selected.columns[2:]),
        "genes": expected,
        "gene_index": {g: i for i, g in enumerate(expected)},
        "x": x,
    }


def load_scanb_raw(universe_set: set[str], primary_titles: list[str]) -> dict:
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
        raise RuntimeError("No frozen-universe genes found in SCAN-B.")

    df = pd.concat(parts, ignore_index=True)
    df = df.groupby("Hugo_Symbol", sort=False, as_index=True).mean(numeric_only=True)

    x = df.to_numpy(dtype=np.float64).T
    finite = np.isfinite(x).all(axis=0)
    xmin = np.full(x.shape[1], np.nan)
    xmax = np.full(x.shape[1], np.nan)
    xmin[finite] = np.min(x[:, finite], axis=0)
    xmax[finite] = np.max(x[:, finite], axis=0)
    eligible = finite & (xmax > xmin)

    genes_all = list(df.index)
    genes = [g for g, ok in zip(genes_all, eligible) if ok]
    x = x[:, eligible]

    if len(genes) != 9220:
        raise RuntimeError(f"Expected 9,220 corrected SCAN-B genes, found {len(genes)}.")

    return {
        "samples": primary_titles,
        "genes": genes,
        "gene_index": {g: i for i, g in enumerate(genes)},
        "x": x,
    }


def source_loading_signs(
    weights: pd.DataFrame,
    program_id: str,
    genes: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    w = weights.loc[weights["program_id"] == program_id].copy()
    w["Hugo_Symbol"] = w["Hugo_Symbol"].map(canon_symbol)
    w = w.set_index("Hugo_Symbol", drop=False)

    loading = np.array(
        [float(w.loc[g, "source_pc1_loading"]) for g in genes],
        dtype=np.float64,
    )
    sign = np.array(
        [float(w.loc[g, "source_pc1_loading_sign"]) for g in genes],
        dtype=np.float64,
    )

    if not np.isfinite(loading).all() or not np.isfinite(sign).all():
        raise RuntimeError(f"{program_id}: non-finite frozen source loading/sign.")

    return loading, sign


def fixed_edge_subset(program_id: str) -> tuple[np.ndarray, np.ndarray]:
    path = (
        SPECIFICITY_ROOT
        / "per_program"
        / "SCANB_GSE96058"
        / f"{program_id}_mapping_specificity_null_v2.npz"
    )
    require(path)
    z = np.load(path)

    return (
        np.asarray(z["edge_slot_i_0based"], dtype=np.int64),
        np.asarray(z["edge_slot_j_0based"], dtype=np.int64),
    )


def platform_bootstrap_rng(
    platform: str,
    program_id: str,
    attempt_id: int,
) -> np.random.Generator:
    ss = np.random.SeedSequence(
        [
            BASE_SEED,
            PLATFORM_CODE[platform],
            module_number(program_id),
            101,
            int(attempt_id),
        ]
    )
    return np.random.default_rng(ss)


def pc1_random_start(
    platform: str,
    program_id: str,
    m: int,
) -> np.ndarray:
    ss = np.random.SeedSequence(
        [
            BASE_SEED,
            PLATFORM_CODE[platform],
            module_number(program_id),
            301,
        ]
    )
    return np.random.default_rng(ss).normal(size=m)


def support_variable_torch(
    x_t: torch.Tensor,
    counts_t: torch.Tensor,
) -> tuple[bool, np.ndarray]:
    support = torch.nonzero(counts_t > 0, as_tuple=False).flatten()
    xs = x_t.index_select(0, support)
    xmin = torch.amin(xs, dim=0)
    xmax = torch.amax(xs, dim=0)
    bad_t = (
        (~torch.isfinite(xmin))
        | (~torch.isfinite(xmax))
        | (xmax <= xmin)
    )
    bad = bad_t.detach().cpu().numpy().astype(bool)
    ok = not bool(torch.any(bad_t))

    del support, xs, xmin, xmax, bad_t
    return ok, bad


def selected_edge_corr_weighted(
    x_t: torch.Tensor,
    counts_t: torch.Tensor,
    ii_t: torch.Tensor,
    jj_t: torch.Tensor,
) -> tuple[torch.Tensor | None, np.ndarray]:
    ok, bad = support_variable_torch(x_t, counts_t)
    if not ok:
        return None, bad

    total = torch.sum(counts_t)
    w = counts_t / total
    mu = torch.sum(w[:, None] * x_t, dim=0)
    xc = x_t - mu[None, :]
    var = torch.sum(w[:, None] * xc * xc, dim=0)

    if bool(torch.any(~torch.isfinite(var))) or bool(torch.any(var <= 0)):
        bad2 = (
            (~torch.isfinite(var))
            | (var <= 0)
        ).detach().cpu().numpy().astype(bool)
        return None, bad2

    e = ii_t.numel()
    out = torch.empty(e, dtype=torch.float64, device=x_t.device)

    for start in range(0, e, EDGE_BATCH):
        stop = min(e, start + EDGE_BATCH)
        ib = ii_t[start:stop]
        jb = jj_t[start:stop]

        xi = xc.index_select(1, ib)
        xj = xc.index_select(1, jb)
        cov = torch.sum(w[:, None] * xi * xj, dim=0)
        out[start:stop] = cov / torch.sqrt(var[ib] * var[jb])

        del ib, jb, xi, xj, cov

    del w, mu, xc, var
    return out, bad


def explicit_edges(
    x: np.ndarray,
    idx: np.ndarray,
    ii: np.ndarray,
    jj: np.ndarray,
    k: int,
) -> np.ndarray:
    xb = x[idx]
    out = np.empty(k, dtype=np.float64)

    for q in range(k):
        out[q] = np.corrcoef(
            xb[:, ii[q]],
            xb[:, jj[q]],
        )[0, 1]

    return out


def cell_paths(platform: str, program_id: str) -> dict[str, Path]:
    d = OUT_DIR / "per_platform" / platform
    d.mkdir(parents=True, exist_ok=True)

    return {
        "checkpoint": d / f"{program_id}_edge_bootstrap_checkpoint_v1.npz",
        "result": d / f"{program_id}_platform_sensitivity_result_v1.json",
        "degenerate": d / f"{program_id}_platform_degenerate_genes_v1.tsv",
    }


def run_platform_cell(
    platform: str,
    program_id: str,
    source_x: np.ndarray,
    target_x: np.ndarray,
    genes: list[str],
    source_loading: np.ndarray,
    source_sign: np.ndarray,
    ii: np.ndarray,
    jj: np.ndarray,
    pooled_row: dict,
    device: torch.device,
) -> dict:
    paths = cell_paths(platform, program_id)

    if paths["result"].exists():
        try:
            old = json.loads(paths["result"].read_text(encoding="utf-8"))
        except Exception:
            old = {}
        if (
            old.get("script_version") == SCRIPT_VERSION
            and old.get("status") in {
                "COMPLETE",
                "NOT_ESTIMABLE_EXACT_VARIATION",
                "BOOTSTRAP_CI_NOT_ESTIMABLE",
            }
        ):
            print("      completed checkpoint found — reusing")
            return old

    bad_target = ~exact_variable_np(target_x)

    if bad_target.any():
        pd.DataFrame(
            {
                "Hugo_Symbol": genes,
                "target_exact_variable": (~bad_target).astype(int),
            }
        ).loc[
            lambda d: d["target_exact_variable"] == 0
        ].to_csv(
            paths["degenerate"],
            sep="\t",
            index=False,
        )

        result = {
            "script_version": SCRIPT_VERSION,
            "status": "NOT_ESTIMABLE_EXACT_VARIATION",
            "platform": platform,
            "program_id": program_id,
            "source_n": int(source_x.shape[0]),
            "target_n": int(target_x.shape[0]),
            "n_fixed_genes": int(len(genes)),
            "target_exact_constant_genes": int(bad_target.sum()),
            "rho_edge_all_edges": None,
            "rho_load_descriptive": None,
            "formal_p_value": None,
            "primary_classification_changed": False,
        }
        paths["result"].write_text(
            json.dumps(result, indent=2),
            encoding="utf-8",
        )
        return result

    source_corr, source_diff = full_corr_from_x(source_x, device)
    target_corr, target_diff = full_corr_from_x(target_x, device)

    tri_i, tri_j = np.triu_indices(len(genes), k=1)
    source_edges = source_corr[tri_i, tri_j]
    target_edges = target_corr[tri_i, tri_j]
    rho_edge = spearman_cpu(source_edges, target_edges)

    target_loading, pcmeta = pc1_loading(
        corr=target_corr,
        x=target_x,
        source_loading=source_loading,
        source_sign=source_sign,
        random_start=pc1_random_start(platform, program_id, len(genes)),
        device=device,
    )
    rho_load = spearman_cpu(source_loading, target_loading)

    source_selected = source_corr[ii, jj]
    source_selected_rank, source_selected_norm = centered_rank_cpu(source_selected)

    del source_corr, target_corr, source_edges, target_edges, target_loading
    gc.collect()

    target_t = torch.as_tensor(
        target_x,
        dtype=torch.float64,
        device=device,
    )
    ii_t = torch.as_tensor(
        ii,
        dtype=torch.long,
        device=device,
    )
    jj_t = torch.as_tensor(
        jj,
        dtype=torch.long,
        device=device,
    )
    source_rank_t = torch.as_tensor(
        source_selected_rank,
        dtype=torch.float64,
        device=device,
    )
    source_norm_t = torch.as_tensor(
        source_selected_norm,
        dtype=torch.float64,
        device=device,
    )

    rho_boot = np.full(
        BOOTSTRAP_REPLICATES,
        np.nan,
        dtype=np.float64,
    )
    completed = 0
    attempts = 0
    degenerate_counts = np.zeros(
        len(genes),
        dtype=np.int64,
    )

    if paths["checkpoint"].exists():
        ck = np.load(paths["checkpoint"])
        completed = int(ck["completed"])
        attempts = int(ck["attempts"])
        rho_boot[:] = ck["rho_boot"]
        degenerate_counts[:] = ck["degenerate_counts"]
        print(
            f"      resuming edge bootstrap valid={completed}/{BOOTSTRAP_REPLICATES}; "
            f"attempts={attempts}"
        )

    engine_validated = completed > 0
    t0 = time.perf_counter()

    while completed < BOOTSTRAP_REPLICATES and attempts < MAX_ATTEMPTS:
        attempts += 1
        rng = platform_bootstrap_rng(
            platform,
            program_id,
            attempts,
        )
        idx = rng.integers(
            0,
            target_x.shape[0],
            size=target_x.shape[0],
            endpoint=False,
        )
        counts = np.bincount(
            idx,
            minlength=target_x.shape[0],
        ).astype(np.float64)
        counts_t = torch.as_tensor(
            counts,
            dtype=torch.float64,
            device=device,
        )

        target_edges_t, bad = selected_edge_corr_weighted(
            target_t,
            counts_t,
            ii_t,
            jj_t,
        )

        if target_edges_t is None:
            degenerate_counts += bad.astype(np.int64)
            del counts_t
            continue

        if not engine_validated:
            k = min(CPU_VALIDATION_EDGES, len(ii))
            ref = explicit_edges(
                target_x,
                idx,
                ii,
                jj,
                k,
            )
            got = target_edges_t[:k].detach().cpu().numpy()
            diff = float(np.max(np.abs(ref - got)))

            if diff > CPU_GPU_TOL:
                raise RuntimeError(
                    f"{platform} {program_id}: bootstrap selected-edge "
                    f"CPU/GPU mismatch {diff:.3e}"
                )

            print(
                f"      bootstrap selected-edge CPU↔GPU validation "
                f"max|Δ|={diff:.3e}"
            )
            engine_validated = True

        rho_boot[completed] = gpu_spearman_against_fixed(
            source_rank_t,
            source_norm_t,
            target_edges_t,
        )
        completed += 1

        del counts_t, target_edges_t

        if (
            completed % CHECKPOINT_EVERY == 0
            or completed == BOOTSTRAP_REPLICATES
        ):
            np.savez_compressed(
                paths["checkpoint"],
                completed=np.array(completed, dtype=np.int64),
                attempts=np.array(attempts, dtype=np.int64),
                rho_boot=rho_boot,
                degenerate_counts=degenerate_counts,
            )
            print(
                f"      edge bootstrap {completed:4d}/{BOOTSTRAP_REPLICATES}; "
                f"attempts={attempts:5d}; "
                f"median rho={np.nanmedian(rho_boot[:completed]):+.4f}; "
                f"elapsed={time.perf_counter()-t0:.1f}s"
            )

    deg_rows = [
        {
            "Hugo_Symbol": gene,
            "degenerate_bootstrap_attempts": int(count),
        }
        for gene, count in zip(genes, degenerate_counts)
        if count > 0
    ]
    pd.DataFrame(
        deg_rows,
        columns=[
            "Hugo_Symbol",
            "degenerate_bootstrap_attempts",
        ],
    ).to_csv(
        paths["degenerate"],
        sep="\t",
        index=False,
    )

    ci_ok = completed == BOOTSTRAP_REPLICATES

    result = {
        "script_version": SCRIPT_VERSION,
        "status": "COMPLETE" if ci_ok else "BOOTSTRAP_CI_NOT_ESTIMABLE",
        "platform": platform,
        "program_id": program_id,
        "source_n": int(source_x.shape[0]),
        "target_n": int(target_x.shape[0]),
        "n_fixed_genes": int(len(genes)),
        "n_all_edges": int(len(tri_i)),
        "n_bootstrap_edges": int(len(ii)),
        "rho_edge_all_edges": rho_edge,
        "rho_load_descriptive": rho_load,
        "pooled_scanb_rho_edge_reference": float(pooled_row["rho_edge"]),
        "pooled_scanb_rho_load_reference": float(pooled_row["rho_load"]),
        "source_full_corr_cpu_gpu_max_abs_diff": source_diff,
        "target_full_corr_cpu_gpu_max_abs_diff": target_diff,
        "target_pc1": pcmeta,
        "edge_bootstrap_valid_replicates": int(completed),
        "edge_bootstrap_attempts": int(attempts),
        "edge_bootstrap_invalid_attempts": int(attempts - completed),
        "edge_bootstrap_median": (
            float(np.median(rho_boot[:completed]))
            if completed
            else None
        ),
        "edge_bootstrap_ci_low": (
            float(np.quantile(rho_boot[:completed], 0.025))
            if completed
            else None
        ),
        "edge_bootstrap_ci_high": (
            float(np.quantile(rho_boot[:completed], 0.975))
            if completed
            else None
        ),
        "formal_p_value": None,
        "primary_classification_changed": False,
    }

    paths["result"].write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )

    del (
        target_t,
        ii_t,
        jj_t,
        source_rank_t,
        source_norm_t,
    )
    torch.cuda.empty_cache()

    return result


def main() -> None:
    print("=" * 150)
    print("Paper 4 / TCBB - SCAN-B platform-specific biological transport sensitivity")
    print("=" * 150)
    print(f"Script version: {SCRIPT_VERSION}")
    print()
    print("Frozen scientific contract:")
    print(f"  Platform groups:                                   {PLATFORM_COUNTS}")
    print("  Source:                                            full frozen 1,082-sample TCGA")
    print("  Program genes:                                     exact corrected 05a v2 SCAN-B sets")
    print("  rho_edge point estimate:                           ALL fixed-gene edges")
    print("  rho_load:                                          descriptive point estimate")
    print(f"  Target edge bootstraps/platform/program:           {BOOTSTRAP_REPLICATES:,}")
    print("  Bootstrap edge subset:                             exact corrected 05b v2 subset")
    print("  Platform-specific gene dropping:                   NO")
    print("  New discovery p-values / classifier changes:       NO")
    print("=" * 150)

    for p in [
        EXEC_CONTRACT,
        PLATFORM_MANIFEST,
        PLATFORM_COMPOSITION,
        CLASSIFICATION,
        UNIVERSE,
        WEIGHTS,
        TCGA_EXPR,
        SCANB_EXPR,
        SCANB_PRIMARY,
    ]:
        require(p)

    ex = json.loads(EXEC_CONTRACT.read_text(encoding="utf-8"))
    if ex.get("status") != "FROZEN_BEFORE_SCANB_PLATFORM_SENSITIVITY_RESULTS":
        raise RuntimeError("05h1 execution contract has unexpected status.")

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable.")
    device = torch.device("cuda:0")
    torch.set_default_dtype(torch.float64)

    prop = torch.cuda.get_device_properties(0)
    print(
        f"CUDA device: {prop.name}; VRAM={prop.total_memory/1024**3:.2f} GB; "
        f"capability={prop.major}.{prop.minor}"
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    universe = pd.read_csv(
        UNIVERSE,
        sep="\t",
        dtype=str,
    ).fillna("")
    if len(universe) != 10000:
        raise RuntimeError(f"Frozen source universe is not 10,000 genes: {len(universe)}")
    universe_set = {
        canon_symbol(x)
        for x in universe["Hugo_Symbol"]
    }

    weights = pd.read_csv(
        WEIGHTS,
        sep="\t",
        dtype=str,
    ).fillna("")

    primary = pd.read_csv(
        SCANB_PRIMARY,
        sep="\t",
        dtype=str,
    ).fillna("")
    primary_titles = primary["primary_title"].tolist()
    if len(primary_titles) != 3273:
        raise RuntimeError("SCAN-B primary manifest is not 3,273 profiles.")

    platform_manifest = pd.read_csv(
        PLATFORM_MANIFEST,
        sep="\t",
        dtype=str,
    ).fillna("")
    if len(platform_manifest) != 3273:
        raise RuntimeError("SCAN-B platform manifest is not 3,273 profiles.")

    if set(platform_manifest["sample_title"]) != set(primary_titles):
        raise RuntimeError("SCAN-B platform manifest does not match frozen primary titles.")

    got_counts = platform_manifest["platform_id"].value_counts().to_dict()
    got_counts = {str(k): int(v) for k, v in got_counts.items()}
    if got_counts != PLATFORM_COUNTS:
        raise RuntimeError(
            f"Platform-count replay mismatch: {got_counts} vs {PLATFORM_COUNTS}"
        )

    print("\n[1/4] Loading frozen source and corrected SCAN-B matrices ...")
    source = load_source_raw(universe)
    scanb = load_scanb_raw(universe_set, primary_titles)
    print(
        f"  TCGA:   {source['x'].shape[0]:,} x {source['x'].shape[1]:,}"
    )
    print(
        f"  SCAN-B: {scanb['x'].shape[0]:,} x {scanb['x'].shape[1]:,}"
    )

    classification = pd.read_csv(
        CLASSIFICATION,
        sep="\t",
        low_memory=False,
    )
    scanb_cls = classification.loc[
        classification["target"] == "SCANB_GSE96058"
    ].copy()

    platform_indices = {}
    sample_index = {
        s: i
        for i, s in enumerate(scanb["samples"])
    }

    for platform, expected_n in PLATFORM_COUNTS.items():
        titles = platform_manifest.loc[
            platform_manifest["platform_id"] == platform,
            "sample_title",
        ].tolist()

        idx = np.array(
            [sample_index[t] for t in titles],
            dtype=np.int64,
        )
        if len(idx) != expected_n:
            raise RuntimeError(
                f"{platform}: expected {expected_n} profiles, found {len(idx)}"
            )
        platform_indices[platform] = idx

    all_results = []
    degenerate_summary = []

    print("\n[2/4] Platform-specific transport effects + edge bootstrap ...")

    assessable_rows = [
        r
        for r in scanb_cls.to_dict(orient="records")
        if boolish(r["primary_assessable"])
    ]

    for pidx, row in enumerate(assessable_rows, start=1):
        program_id = str(row["program_id"])

        gene_file = (
            DIRECT_ROOT
            / "per_program"
            / "SCANB_GSE96058"
            / f"{program_id}_evaluable_genes_v2.tsv"
        )
        require(gene_file)
        gdf = pd.read_csv(
            gene_file,
            sep="\t",
            dtype=str,
        ).fillna("")
        genes = [
            canon_symbol(x)
            for x in gdf["Hugo_Symbol"]
        ]

        source_cols = np.array(
            [source["gene_index"][g] for g in genes],
            dtype=np.int64,
        )
        target_cols = np.array(
            [scanb["gene_index"][g] for g in genes],
            dtype=np.int64,
        )

        source_x = source["x"][:, source_cols]
        source_loading, source_sign = source_loading_signs(
            weights,
            program_id,
            genes,
        )

        ii, jj = fixed_edge_subset(program_id)
        if int(max(ii.max(), jj.max())) >= len(genes):
            raise RuntimeError(
                f"{program_id}: frozen edge subset exceeds corrected gene set."
            )

        print(
            f"\n  [{pidx:02d}/{len(assessable_rows):02d}] {program_id}: "
            f"genes={len(genes):,}; bootstrap edges={len(ii):,}"
        )

        for platform in ["GPL11154", "GPL18573"]:
            idx = platform_indices[platform]
            target_x = scanb["x"][idx[:, None], target_cols[None, :]]

            print(
                f"    {platform}: target n={target_x.shape[0]:,}"
            )

            result = run_platform_cell(
                platform=platform,
                program_id=program_id,
                source_x=source_x,
                target_x=target_x,
                genes=genes,
                source_loading=source_loading,
                source_sign=source_sign,
                ii=ii,
                jj=jj,
                pooled_row=row,
                device=device,
            )
            all_results.append(result)

            if result["status"] == "COMPLETE":
                print(
                    f"      rho_edge(all)={result['rho_edge_all_edges']:+.4f}; "
                    f"bootstrap median={result['edge_bootstrap_median']:+.4f} "
                    f"[{result['edge_bootstrap_ci_low']:+.4f},"
                    f"{result['edge_bootstrap_ci_high']:+.4f}]; "
                    f"rho_load={result['rho_load_descriptive']:+.4f}; "
                    f"attempts={result['edge_bootstrap_attempts']}"
                )
            else:
                print(
                    f"      STATUS={result['status']}"
                )

    print("\n[3/4] Replaying frozen platform × PAM50 composition ...")

    comp = pd.read_csv(
        PLATFORM_COMPOSITION,
        sep="\t",
        low_memory=False,
    )
    comp_counts = (
        comp.groupby("platform_id")["n"]
        .sum()
        .astype(int)
        .to_dict()
    )
    for platform, expected_n in PLATFORM_COUNTS.items():
        if int(comp_counts.get(platform, -1)) != expected_n:
            raise RuntimeError(
                f"{platform}: PAM50 composition does not sum to {expected_n}."
            )
    print("  platform × PAM50 composition replay: PASS")

    print("\n[4/4] Writing platform-sensitivity outputs ...")

    result_df = pd.DataFrame(all_results)
    result_path = OUT_DIR / "scanb_platform_sensitivity_summary_v1.tsv"
    result_df.to_csv(
        result_path,
        sep="\t",
        index=False,
    )

    comp_out = OUT_DIR / "scanb_platform_pam50_composition_replay_v1.tsv"
    comp.to_csv(
        comp_out,
        sep="\t",
        index=False,
    )

    complete = int((result_df["status"] == "COMPLETE").sum())
    noncomplete = int(len(result_df) - complete)

    per_platform = {}
    for platform in ["GPL11154", "GPL18573"]:
        d = result_df.loc[result_df["platform"] == platform]
        per_platform[platform] = {
            "programs": int(len(d)),
            "complete": int((d["status"] == "COMPLETE").sum()),
            "not_complete": int((d["status"] != "COMPLETE").sum()),
        }

    master = {
        "script_version": SCRIPT_VERSION,
        "status": "SCANB_PLATFORM_SENSITIVITY_COMPLETE",
        "platform_counts": PLATFORM_COUNTS,
        "primary_assessable_programs": int(len(assessable_rows)),
        "platform_program_cells": int(len(result_df)),
        "complete_cells": complete,
        "noncomplete_cells": noncomplete,
        "per_platform": per_platform,
        "edge_bootstrap_valid_replicates_required": BOOTSTRAP_REPLICATES,
        "loading_bootstrap_performed": False,
        "formal_new_discovery_p_values": False,
        "primary_classification_changed": False,
        "target_outcomes_or_treatment_loaded": False,
        "summary_file": str(result_path),
        "pam50_composition_file": str(comp_out),
    }

    master_path = OUT_DIR / "scanb_platform_sensitivity_v1.json"
    master_path.write_text(
        json.dumps(master, indent=2),
        encoding="utf-8",
    )

    print("\n" + "=" * 150)
    print("05h SCAN-B PLATFORM SENSITIVITY: COMPLETE")
    print("=" * 150)
    print(
        f"Platform/program COMPLETE cells:       "
        f"{complete}/{len(result_df)}"
    )
    print(
        f"GPL11154 COMPLETE:                     "
        f"{per_platform['GPL11154']['complete']}/"
        f"{per_platform['GPL11154']['programs']}"
    )
    print(
        f"GPL18573 COMPLETE:                     "
        f"{per_platform['GPL18573']['complete']}/"
        f"{per_platform['GPL18573']['programs']}"
    )
    print("PAM50 composition replay:              PASS")
    print("New discovery p-values:                NO")
    print("Primary classifications changed:       NO")
    print()
    print(f"Summary:     {result_path}")
    print(f"Composition: {comp_out}")
    print(f"Master:      {master_path}")
    print("=" * 150)

    del source, scanb
    torch.cuda.empty_cache()
    gc.collect()


if __name__ == "__main__":
    main()
